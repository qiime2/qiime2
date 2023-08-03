import os
import pathlib
import shutil
import tempfile
import unittest
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Generator

from ..parse import ProvDAG
from ..util import UUID

from qiime2 import Artifact
from qiime2.sdk.plugin_manager import PluginManager


@dataclass
class TestArtifact:
    name: str
    artifact: Artifact
    uuid: str
    filepath: str
    dag: ProvDAG
    archive_version: int = 6


class TestArtifacts:
    def __init__(self):
        self.pm = PluginManager()
        self.dp = self.pm.plugins['dummy-plugin']
        self.tempdir = tempfile.mkdtemp(prefix='qiime2-dummy-artifacts-temp-')

        # TODO: move versioned artifacts into root of data dir once everything
        # else is gone
        self.datadir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'data',
            'versioned-artifacts'
        )

        # artifacts with no action provenance (except import)
        single_int = Artifact.import_data('SingleInt', 0)
        int_seq1 = Artifact.import_data('IntSequence1', [1, 1, 2])
        int_seq2 = Artifact.import_data('IntSequence2', [3, 5])
        mapping1 = Artifact.import_data('Mapping', {'a': 42})
        mapping2 = Artifact.import_data('Mapping', {'c': 8, 'd': 13})

        for artifact, name in zip(
            [single_int, int_seq1, int_seq2, mapping1, mapping2],
            ['single_int', 'int_seq1', 'int_seq2', 'mapping1', 'mapping2']
        ):
            fp = os.path.join(self.tempdir, f'{name}.qza')
            artifact.save(fp)
            test_artifact = TestArtifact(
                name, artifact, str(artifact.uuid), fp, ProvDAG(fp)
            )
            setattr(self, name, test_artifact)

        concat_ints = self.dp.methods['concatenate_ints']
        split_ints = self.dp.methods['split_ints']
        merge_mappings = self.dp.methods['merge_mappings']
        # to be used
        # identity_with_metadata = self.dp.actions['identity_with_metadata']

        # artifacts with some simple actions in their provenance
        concated_ints, = concat_ints(int_seq1, int_seq1, int_seq2, 7, 13)
        other_concated_ints, = concat_ints(int_seq1, int_seq1, int_seq2,
                                           81, 64)
        splitted_ints, _ = split_ints(int_seq2)
        merged_mappings, = merge_mappings(mapping1, mapping2)

        # create artifact with pipeline provenance
        typical_pipeline = self.dp.pipelines['typical_pipeline']
        _, _, _, pipeline_viz, _ = typical_pipeline(int_seq1, mapping1, False)

        for artifact, name in zip(
            [concated_ints, other_concated_ints, splitted_ints,
             merged_mappings, pipeline_viz],
            ['concated_ints', 'other_concated_ints', 'splitted_ints',
             'merged_mappings', 'pipeline_viz']
        ):
            if name == 'pipeline_viz':
                ext = '.qzv'
            else:
                ext = '.qza'

            fp = os.path.join(self.tempdir, f'{name}{ext}')
            artifact.save(fp)
            test_artifact = TestArtifact(
                name, artifact, str(artifact.uuid), fp, ProvDAG(fp)
            )
            setattr(self, name, test_artifact)

        # import dummy artifacts for versions [0, 6]
        for version in range(0, 7):
            if version == 0:
                dirname = 'table-v0'
            else:
                dirname = f'concated-ints-v{version}'

            versioned_artifact_dir = os.path.join(self.datadir, dirname)
            temp_zf_path = os.path.join(self.tempdir, 'temp.zip')
            zf = zipfile.ZipFile(temp_zf_path, 'w', zipfile.ZIP_DEFLATED)
            for root, dirs, files in os.walk(versioned_artifact_dir):
                for file in files:
                    filepath = os.path.join(root, file)
                    zf.write(
                        filepath,
                        os.path.relpath(filepath, versioned_artifact_dir)
                    )
            zf.close()

            filename = f'{dirname}.qza'
            fp = os.path.join(self.tempdir, filename)

            if version == 0:
                shutil.copy(temp_zf_path, fp)
                a = None
            else:
                a = Artifact.load(temp_zf_path)
                a.save(fp)

            dag = ProvDAG(fp)
            assert len(dag.terminal_nodes) == 1
            terminal_node, *_ = dag.terminal_nodes
            uuid = terminal_node._uuid

            name = filename.replace('-', '_').replace('.qza', '')
            ta = TestArtifact(name, a, uuid, fp, dag, version)
            setattr(self, name, ta)

        # create archive with missing checksums.md5
        with generate_archive_with_file_removed(
            self.single_int.filepath,
            self.single_int.uuid,
            'checksums.md5'
        ) as altered_archive:
            self.dag_missing_md5 = ProvDAG(altered_archive)

    def free(self):
        shutil.rmtree(self.tempdir)


class CustomAssertions(unittest.TestCase):
    def assertREAppearsOnlyOnce(self, text, only_once, msg=None):
        appears_once_re = \
            (f'(?s)^(?:(?!{only_once}).)*{only_once}(?!.*{only_once}).*$')
        self.assertRegex(text, appears_once_re, msg)


def is_root_provnode_data(fp):
    """
    a filter predicate which returns metadata, action, citation,
    and VERSION fps with which we can construct a ProvNode
    """
    # Handle provenance files...
    if 'provenance' in fp and 'artifacts' not in fp \
        and ('action.yaml' in fp or
             'citations.bib' in fp
             ):
        return True

    # then handle files available at root, which require a cast
    if pathlib.Path(fp).parts[1] in ('VERSION',
                                     'metadata.yaml',
                                     'checksums.md5'):
        return True


@contextmanager
def generate_archive_with_file_removed(qzv_fp: str, root_uuid: UUID,
                                       file_to_drop: pathlib.Path) -> \
                                           Generator[pathlib.Path, None, None]:
    """
    Deleting files from zip archives is hard, so this makes a temporary
    copy of qzf_fp with fp_to_drop removed and returns a handle to this archive

    file_to_drop should represent the relative path to the file within the
    zip archive, excluding the root directory (named for the root UUID).

    e.g. `/d9e080bb-e245-4ab0-a2cf-0a89b63b8050/metadata.yaml` should be passed
    in as `metadata.yaml`

    adapted from https://stackoverflow.com/a/513889/9872253
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_arc = pathlib.Path(tmpdir) / 'mangled.qzv'
        fp_pfx = pathlib.Path(root_uuid)
        zin = zipfile.ZipFile(qzv_fp, 'r')
        zout = zipfile.ZipFile(str(tmp_arc), 'w')
        for item in zin.infolist():
            buffer = zin.read(item.filename)
            drop_filename = str(fp_pfx / file_to_drop)
            if (item.filename != drop_filename):
                zout.writestr(item, buffer)
        zout.close()
        zin.close()
        yield tmp_arc


class ReallyEqualMixin(object):
    """
    Mixin for testing implementations of __eq__/__ne__.

    Based on this public domain code (also explains why the mixin is useful):
    https://ludios.org/testing-your-eq-ne-cmp/
    """

    def assertReallyEqual(self, a, b):
        # assertEqual first, because it will have a good message if the
        # assertion fails.
        self.assertEqual(a, b)
        self.assertEqual(b, a)
        self.assertTrue(a == b)
        self.assertTrue(b == a)
        self.assertFalse(a != b)
        self.assertFalse(b != a)

    def assertReallyNotEqual(self, a, b):
        # assertNotEqual first, because it will have a good message if the
        # assertion fails.
        self.assertNotEqual(a, b)
        self.assertNotEqual(b, a)
        self.assertFalse(a == b)
        self.assertFalse(b == a)
        self.assertTrue(a != b)
        self.assertTrue(b != a)
