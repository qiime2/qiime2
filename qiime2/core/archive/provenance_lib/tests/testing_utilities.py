import os
import pathlib
import shutil
import tempfile
import unittest
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Generator
from zipfile import ZipFile, ZIP_DEFLATED

from ..parse import ProvDAG

import qiime2
from qiime2 import Artifact, Metadata, ResultCollection
from qiime2.core.archive import Archiver
from qiime2.sdk.plugin_manager import PluginManager
from qiime2.plugin import Plugin, Int
from qiime2.core.testing.type import IntSequence1
from qiime2.core.testing.method import concatenate_ints


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

        self.datadir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'data'
        )

        self.init_import_artifacts()
        self.init_action_artifacts()
        self.init_all_version_artifacts()
        self.init_artifact_with_md_in_provenance()
        self.init_no_checksum_dag()

    def init_import_artifacts(self):
        '''
        artifacts with only import in their provenance
        '''
        single_int = Artifact.import_data('SingleInt', 0)
        single_int2 = Artifact.import_data('SingleInt', 7)
        int_seq1 = Artifact.import_data('IntSequence1', [1, 1, 2])
        int_seq2 = Artifact.import_data('IntSequence2', [3, 5])
        mapping1 = Artifact.import_data('Mapping', {'a': 42})
        mapping2 = Artifact.import_data('Mapping', {'c': 8, 'd': 13})

        for name in (
            'single_int', 'single_int2', 'int_seq1', 'int_seq2', 'mapping1',
            'mapping2'
        ):
            artifact = locals()[name]
            fp = os.path.join(self.tempdir, f'{name}.qza')
            artifact.save(fp)
            test_artifact = TestArtifact(
                name, artifact, str(artifact.uuid), fp, ProvDAG(fp)
            )
            setattr(self, name, test_artifact)

    def init_action_artifacts(self):
        '''
        artifacts that have at least one non-import action in their provenance
        '''
        concat_ints = self.dp.methods['concatenate_ints']
        split_ints = self.dp.methods['split_ints']
        merge_mappings = self.dp.methods['merge_mappings']
        identity_with_metadata = self.dp.actions['identity_with_metadata']
        identity_with_metadata_column = \
            self.dp.actions['identity_with_metadata_column']
        dict_of_ints = self.dp.actions['dict_of_ints']
        optional_artifacts_method = \
            self.dp.actions['optional_artifacts_method']

        concated_ints, = concat_ints(
            self.int_seq1.artifact, self.int_seq1.artifact,
            self.int_seq2.artifact, 7, 13
        )
        other_concated_ints, = concat_ints(
            self.int_seq1.artifact, self.int_seq1.artifact,
            self.int_seq2.artifact, 81, 64
        )
        splitted_ints, _ = split_ints(self.int_seq2.artifact)
        merged_mappings, = merge_mappings(
            self.mapping1.artifact, self.mapping2.artifact
        )

        # artifacts with input artifact viewed as metadata
        int_seq_with_md, = identity_with_metadata(
            self.int_seq1.artifact, self.mapping1.artifact.view(Metadata)
        )
        int_seq_with_md_column, = identity_with_metadata_column(
            self.int_seq1.artifact,
            self.mapping2.artifact.view(Metadata).get_column('c')
        )
        concated_ints_with_md_column, = concat_ints(
            self.int_seq1.artifact, int_seq_with_md_column,
            self.int_seq2.artifact, 69, 2001
        )

        # artifact with input collection
        ints_dict = {
            'int1': self.single_int.artifact,
            'int2': self.single_int2.artifact,
        }
        ints_collection = ResultCollection(ints_dict)
        ints_from_collection, = dict_of_ints(ints_collection)
        int_from_collection = ints_from_collection['int1']

        # artifact with optional inputs left to default None
        int_seq_optional_input, = optional_artifacts_method(
            self.int_seq1.artifact, 8
        )

        # artifact from pipeline
        typical_pipeline = self.dp.pipelines['typical_pipeline']
        _, _, _, pipeline_viz, _ = typical_pipeline(
            self.int_seq1.artifact, self.mapping1.artifact, False
        )

        for name in (
            'concated_ints', 'other_concated_ints', 'splitted_ints',
            'merged_mappings', 'pipeline_viz', 'int_seq_with_md',
            'concated_ints_with_md_column', 'int_from_collection',
            'int_seq_optional_input'
        ):
            artifact = locals()[name]
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

    def init_all_version_artifacts(self):
        '''
        import artifacts for all archive versions (0-6), which are stored
        uncompressed in the test data directory--necessary because we can not
        make artifacts of non-current versions on the fly
        '''
        for version in range(0, 7):
            if version == 0:
                dirname = 'table-v0'
            else:
                dirname = f'concated-ints-v{version}'

            versioned_artifact_dir = os.path.join(self.datadir, dirname)
            temp_zf_path = os.path.join(self.tempdir, 'temp.zip')
            write_zip_file(temp_zf_path, versioned_artifact_dir)

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

    def init_artifact_with_md_in_provenance(self):
        dirname = 'concated-ints-with-md'
        artifact_dir = os.path.join(self.datadir, dirname)
        temp_zf_path = os.path.join(self.tempdir, 'temp.zip')
        write_zip_file(temp_zf_path, artifact_dir)
        filename = f'{dirname}.qza'
        fp = os.path.join(self.tempdir, filename)
        a = Artifact.load(temp_zf_path)
        a.save(fp)

        dag = ProvDAG(fp)
        terminal_node, *_ = dag.terminal_nodes
        uuid = terminal_node._uuid
        name = filename.replace('-', '_').replace('.qza', '')
        ta = TestArtifact(name, a, uuid, fp, dag, 6)
        setattr(self, name, ta)

    def init_no_checksum_dag(self):
        '''
        create archive with missing checksums.md5
        '''
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
            self.concated_ints_v5, self.concated_ints_v6
        )

    def free(self):
        shutil.rmtree(self.tempdir)


def write_zip_file(zfp, unzipped_dir):
    zf = ZipFile(zfp, 'w', ZIP_DEFLATED)
    for root, dirs, files in os.walk(unzipped_dir):
        for file in files:
            filepath = os.path.join(root, file)
            zf.write(
                filepath, os.path.relpath(filepath, unzipped_dir)
            )
    zf.close()


class CustomAssertions(unittest.TestCase):
    def assertREAppearsOnlyOnce(self, text, only_once, msg=None):
        appears_once_re = \
            (f'(?s)^(?:(?!{only_once}).)*{only_once}(?!.*{only_once}).*$')
        self.assertRegex(text, appears_once_re, msg)


def is_root_provnode_data(fp):
    '''
    a filter predicate which returns metadata, action, citation,
    and VERSION fps with which we can construct a ProvNode
    '''
    # Handle provenance files...
    if 'provenance' in fp and 'artifacts' not in fp:
        if 'action.yaml' in fp or 'citations.bib' in fp:
            return True

    # then handle files available at root, which require a cast
    if pathlib.Path(fp).parts[1] in (
        'VERSION', 'metadata.yaml', 'checksums.md5'
    ):
        return True


@contextmanager
def generate_archive_with_file_removed(
    qzv_fp: str, root_uuid: str, file_to_drop: pathlib.Path
) -> Generator[pathlib.Path, None, None]:
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
        zin = ZipFile(qzv_fp, 'r')
        zout = ZipFile(str(tmp_arc), 'w')
        for item in zin.infolist():
            buffer = zin.read(item.filename)
            drop_filename = str(fp_pfx / file_to_drop)
            if (item.filename != drop_filename):
                zout.writestr(item, buffer)
        zout.close()
        zin.close()
        yield tmp_arc


@contextmanager
def monkeypatch_archive_version(patch_version):
    try:
        og_version = Archiver.CURRENT_FORMAT_VERSION
        Archiver.CURRENT_FORMAT_VERSION = patch_version
        yield
    finally:
        Archiver.CURRENT_FORMAT_VERSION = og_version


@contextmanager
def monkeypatch_framework_version(patch_version):
    try:
        og_version = qiime2.__version__
        qiime2.__version__ = patch_version
        yield
    finally:
        qiime2.__version__ = og_version


def write_zip_archive(zfp, unzipped_dir):
    with ZipFile(zfp, 'w') as zf:
        for root, dirs, files in os.walk(unzipped_dir):
            for file in files:
                path = os.path.join(root, file)
                archive_name = os.path.relpath(path, start=unzipped_dir)
                zf.write(path, arcname=archive_name)


other_plugin = Plugin(
    name='other-plugin',
    description='',
    short_description='',
    version='0.0.0-dev',
    website='',
    package='qiime2.core.archive.provenance_lib.tests',
    user_support_text='',
    citations=[]
)
other_plugin.methods.register_function(
    function=concatenate_ints,
    inputs={
        'ints1': IntSequence1,
        'ints2': IntSequence1,
        'ints3': IntSequence1,
    },
    parameters={
        'int1': Int,
        'int2': Int
    },
    outputs={
        'concatenated_ints': IntSequence1
    },
    name='Concatenate integers',
    description='Some description'
)
