import os
import pathlib
import shutil
import tempfile
import unittest
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Generator

from ..archive_parser import (
    ParserV0, ParserV1, ParserV2, ParserV3, ParserV4, ParserV5
)
from .._checksum_validator import ValidationCode
from ..parse import ProvDAG
from ..util import UUID

from qiime2 import Artifact
from qiime2.sdk.plugin_manager import PluginManager
from qiime2.core.archive.archiver import ChecksumDiff


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

    @property
    def all_artifact_versions(self):
        return (
            self.table_v0, self.concated_ints_v1, self.concated_ints_v2,
            self.concated_ints_v3, self.concated_ints_v4,
            self.concated_ints_v5, self.concated_ints_v6,
        )

    def free(self):
        shutil.rmtree(self.tempdir)


DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
TEST_DATA = {
    '0': {'parser': ParserV0,
          'av': '0',
          'fwv': '2.0.5',
          'uuid': '0b8b47bd-f2f8-4029-923c-0e37a68340c3',
          'n_res': 1,
          'qzv_fp': os.path.join(DATA_DIR, 'v0_uu_emperor.qzv'),
          'has_prov': False,
          'prov_is_valid': ValidationCode.PREDATES_CHECKSUMS,
          'checksum': None,
          },
    '1': {'parser': ParserV1,
          'av': '1',
          'fwv': '2.0.6',
          'uuid': '0b8b47bd-f2f8-4029-923c-0e37a68340c3',
          'nonroot_node': '60cde83c-180d-40cb-87c9-b9363f23f796',
          'n_res': 10,
          'qzv_fp': os.path.join(DATA_DIR, 'v1_uu_emperor.qzv'),
          'has_prov': True,
          'prov_is_valid': ValidationCode.PREDATES_CHECKSUMS,
          'checksum': None,
          },
    '2a': {'parser': ParserV2,
           'av': '2',
           'fwv': '2017.9.0',
           'uuid': '219c4bdf-f2b1-4b3f-b66a-08de8a4d17ca',
           'nonroot_node': '512ced83-cc8b-4bed-8c22-a829e8fc89a2',
           'n_res': 10,
           'qzv_fp': os.path.join(DATA_DIR, 'v2a_uu_emperor.qzv'),
           'has_prov': True,
           'prov_is_valid': ValidationCode.PREDATES_CHECKSUMS,
           'checksum': None,
           },
    '2b': {'parser': ParserV2,
           'av': '2',
           'fwv': '2017.10.0',
           'uuid': '8abf8dee-0047-4a7f-9826-e66893182978',
           'nonroot_node': '10ebb316-169e-422c-8fb9-423e131fe42f',
           'n_res': 14,
           'qzv_fp': os.path.join(DATA_DIR, 'v2b_uu_emperor.qzv'),
           'has_prov': True,
           'prov_is_valid': ValidationCode.PREDATES_CHECKSUMS,
           'checksum': None,
           },
    '3': {'parser': ParserV3,
          'av': '3',
          'fwv': '2017.12.0',
          'uuid': '3544061c-6e2f-4328-8345-754416828cb5',
          'nonroot_node': '32c222f5-d991-4168-bca2-d305513e258f',
          'n_res': 14,
          'qzv_fp': os.path.join(DATA_DIR, 'v3_uu_emperor.qzv'),
          'has_prov': True,
          'prov_is_valid': ValidationCode.PREDATES_CHECKSUMS,
          'checksum': None,
          },
    '4': {'parser': ParserV4,
          'av': '4',
          'fwv': '2018.4.0',
          'uuid': '91c2189a-2d2e-4d53-98ee-659caaf6ffc2',
          'nonroot_node': '48c153b4-314c-4249-88a3-020f5444a76f',
          'n_res': 14,
          'qzv_fp': os.path.join(DATA_DIR, 'v4_uu_emperor.qzv'),
          'has_prov': True,
          'prov_is_valid': ValidationCode.PREDATES_CHECKSUMS,
          'checksum': None,
          },
    '5': {'parser': ParserV5,
          'av': '5',
          'fwv': '2018.11.0',
          'uuid': 'ffb7cee3-2f1f-4988-90cc-efd5184ef003',
          'nonroot_node': '3b7d36ff-37ab-4ac2-958b-6a547d442bcf',
          'n_res': 15,
          'qzv_fp': os.path.join(DATA_DIR, 'v5_uu_emperor.qzv'),
          'has_prov': True,
          'prov_is_valid': ValidationCode.VALID,
          'checksum': ChecksumDiff({}, {}, {}),
          },
    }


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
