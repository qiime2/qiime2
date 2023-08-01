from contextlib import redirect_stdout
from dataclasses import dataclass
import io
import networkx as nx
import os
import pathlib
import random
import shutil
import tempfile
import unittest
import warnings
import zipfile

from networkx import DiGraph

from qiime2 import Artifact
from qiime2.sdk.plugin_manager import PluginManager
from qiime2.core.archive.archiver import ChecksumDiff

from .._checksum_validator import ValidationCode
from ..parse import (
    ProvDAG, UnparseableDataError, DirectoryParser, EmptyParser, ProvDAGParser,
    archive_not_parsed, select_parser, parse_provenance,
)
from ..util import UUID
from ..archive_parser import (
    ParserV0, ParserV1, ParserV2, ParserV3, ParserV4, ParserV5,
    Config, ProvNode, ParserResults, ArchiveParser,
)

from .testing_utilities import (
    is_root_provnode_data, generate_archive_with_file_removed,
)

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


@dataclass
class TestArtifact:
    name: str
    artifact: Artifact
    uuid: str
    filepath: str
    dag: ProvDAG
    archive_version: int = 6


class ProvDAGTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dags = dict()
        for archive_version in TEST_DATA:
            # supress warning from parsing provenance for a v0 provDag
            uuid = TEST_DATA['0']['uuid']
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore',  f'Art.*{uuid}.*prior')
                cls.dags[archive_version] = ProvDAG(
                    str(TEST_DATA[archive_version]['qzv_fp']))

    def save_artifact_to_dir(self, artifact, directory):
        dir_fp = os.path.join(self.tempdir, directory)
        try:
            os.mkdir(dir_fp)
        except FileExistsError:
            pass

        fp = os.path.join(dir_fp, f'{artifact.name}.qza')
        artifact.artifact.save(fp)
        return dir_fp

    def setUp(self):
        self.pm = PluginManager()
        self.dp = self.pm.plugins['dummy-plugin']
        self.tempdir = tempfile.mkdtemp(prefix='qiime2-test-parse-temp-')

        # TODO: move versioned artifacts into root of data dir once everything
        # TODO: else is gone
        self.datadir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'data',
            'versioned-artifacts'
        )

        class TestArtifacts:
            # for attaching test artifacts
            pass

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
            setattr(TestArtifacts, name, test_artifact)

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
            setattr(TestArtifacts, name, test_artifact)

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
            setattr(TestArtifacts, name, ta)

        # create archive with missing checksums.md5
        with generate_archive_with_file_removed(
            TestArtifacts.single_int.filepath,
            TestArtifacts.single_int.uuid,
            'checksums.md5'
        ) as altered_archive:
            self.dag_missing_md5 = ProvDAG(altered_archive)

        self.tas = TestArtifacts

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_number_of_nodes(self):
        num_single_int_nodes = len(self.tas.single_int.dag.nodes)
        self.assertEqual(num_single_int_nodes, 1)

        num_int_seq1_nodes = len(self.tas.int_seq1.dag.nodes)
        self.assertEqual(num_int_seq1_nodes, 1)

        num_concated_ints_nodes = len(self.tas.concated_ints.dag.nodes)
        self.assertEqual(num_concated_ints_nodes, 3)

    def test_number_of_nodes_pipeline(self):
        # input int-seq, input mapping, pipeline output viz, aliased viz,
        # split int-seq
        num_pipeline_viz_nodes = len(self.tas.pipeline_viz.dag.nodes)
        self.assertEqual(num_pipeline_viz_nodes, 5)

    def test_number_of_terminal_nodes(self):
        num_int_seq1_term_nodes = len(self.tas.int_seq1.dag.terminal_nodes)
        self.assertEqual(num_int_seq1_term_nodes, 1)

        num_concated_ints_term_nodes = \
            len(self.tas.concated_ints.dag.terminal_nodes)
        self.assertEqual(num_concated_ints_term_nodes, 1)

        self.save_artifact_to_dir(self.tas.int_seq1, 'two-ints')
        dir_fp = self.save_artifact_to_dir(self.tas.int_seq2, 'two-ints')
        dag = ProvDAG(dir_fp)
        self.assertEqual(len(dag.terminal_nodes), 2)

        self.save_artifact_to_dir(self.tas.int_seq1, 'int-and-concated-ints')
        dir_fp = self.save_artifact_to_dir(self.tas.concated_ints,
                                           'int-and-concated-ints')
        dag = ProvDAG(dir_fp)
        self.assertEqual(len(dag.terminal_nodes), 1)

    def test_number_of_terminal_nodes_pipeline(self):
        num_pipeline_viz_term_nodes = \
            len(self.tas.pipeline_viz.dag.terminal_nodes)
        self.assertEqual(num_pipeline_viz_term_nodes, 1)

    def test_root_node_is_archive_root(self):
        with zipfile.ZipFile(self.tas.concated_ints.filepath) as zf:
            all_filenames = zf.namelist()
            root_filenames = filter(is_root_provnode_data, all_filenames)
            root_filepaths = [pathlib.Path(fp) for fp in root_filenames]
            exp_node = ProvNode(Config(), zf, root_filepaths)
            act_terminal_node, *_ = self.tas.concated_ints.dag.terminal_nodes
            self.assertEqual(exp_node, act_terminal_node)

    def test_number_of_actions(self):
        self.assertEqual(self.tas.int_seq1.dag.dag.number_of_edges(), 0)

        self.assertEqual(self.tas.concated_ints.dag.dag.number_of_edges(), 2)

    def test_number_of_actions_pipeline(self):
        # (1) one edge from input intseq to pipeline output viz
        # (2) one edge from input mapping to pipeline output viz
        # (3) one edge from input intseq to left split int
        # (4) one edge from left split int to true output viz, of which
        # pipeline output viz is an alias
        self.assertEqual(self.tas.pipeline_viz.dag.dag.number_of_edges(), 4)

    def test_nonexistent_fp(self):
        fp = os.path.join(self.tempdir, 'does-not-exist.qza')
        err_msg = f'FileNotFoundError.*ArchiveParser.*{fp}'
        with self.assertRaisesRegex(UnparseableDataError, err_msg):
            ProvDAG(fp)

    def test_insufficient_permissions(self):
        fp = self.tas.int_seq1.filepath
        os.chmod(fp, 0o000)
        err_msg = f'PermissionError.*ArchiveParser.*denied.*{fp}'
        with self.assertRaisesRegex(UnparseableDataError, err_msg):
            ProvDAG(fp)

        os.chmod(fp, 0o123)
        err_msg = f'PermissionError.*ArchiveParser.*denied.*{fp}'
        with self.assertRaisesRegex(UnparseableDataError, err_msg):
            ProvDAG(fp)

    def test_not_a_zip_file(self):
        fp = os.path.join(self.tempdir, 'not-a-zip.txt')
        with open(fp, 'w') as fh:
            fh.write("This is just a text file.")

        err_msg = 'zipfile.BadZipFile.*ArchiveParser.*File is not a zip file'
        with self.assertRaisesRegex(UnparseableDataError, err_msg):
            ProvDAG(fp)

    def test_has_digraph(self):
        self.assertIsInstance(self.tas.int_seq1.dag.dag, DiGraph)

        self.assertIsInstance(self.tas.concated_ints.dag.dag, DiGraph)

    def test_dag_attributes(self):
        dag = self.tas.int_seq1.dag
        terminal_node, *_ = dag.terminal_nodes
        self.assertIsInstance(terminal_node, ProvNode)

        self.assertEqual(dag.provenance_is_valid, ValidationCode.VALID)

        empty_checksum_diff = ChecksumDiff(added={}, removed={}, changed={})
        self.assertEqual(dag.checksum_diff, empty_checksum_diff)

    def test_node_action_names(self):
        int_seq1_node, *_ = self.tas.int_seq1.dag.terminal_nodes
        self.assertEqual(int_seq1_node.action.action_name, 'import')

        concated_ints_node, *_ = self.tas.concated_ints.dag.terminal_nodes
        self.assertEqual(concated_ints_node.action.action_name,
                         'concatenate_ints')

    def test_node_action_names_pipeline(self):
        pipeline_viz_node, *_ = self.tas.pipeline_viz.dag.terminal_nodes
        self.assertEqual(pipeline_viz_node.action.action_name,
                         'typical_pipeline')

    def test_has_correct_edges(self):
        edges = self.tas.concated_ints.dag.dag.edges

        self.assertIn(
            (self.tas.int_seq1.uuid, self.tas.concated_ints.uuid), edges
        )
        self.assertIn(
            (self.tas.int_seq2.uuid, self.tas.concated_ints.uuid), edges
        )
        self.assertNotIn(
            (self.tas.int_seq1.uuid, self.tas.int_seq2.uuid), edges
        )
        self.assertNotIn(
            (self.tas.concated_ints.uuid, self.tas.int_seq1.uuid), edges
        )

    def test_dag_repr(self):
        exp_repr = f'ProvDAG.*Artifacts.*{self.tas.int_seq1.uuid}'
        self.assertRegex(repr(self.tas.int_seq1.dag), exp_repr)

    def test_node_repr(self):
        int_seq1_node, *_ = self.tas.int_seq1.dag.terminal_nodes
        uuid = int_seq1_node._uuid
        type_ = int_seq1_node.type
        format_ = int_seq1_node.format

        exp_repr = f'(?s)UUID.*{uuid}.*Type:.*{type_}.*Data Format:.*{format_}'
        self.assertRegex(repr(int_seq1_node), exp_repr)

    def test_dag_eq(self):
        dag = self.tas.int_seq1.dag
        self.assertEqual(dag, dag)

        fp = self.tas.int_seq1.filepath
        self.assertEqual(ProvDAG(fp), ProvDAG(fp))

        # because they are isomorphic
        self.assertEqual(self.tas.int_seq1.dag, self.tas.int_seq2.dag)

    def test_dag_not_eq(self):
        self.assertNotEqual(self.tas.int_seq1.dag, self.tas.concated_ints.dag)

    def test_captures_full_history(self):
        concat_ints = self.dp.actions['concatenate_ints']

        next_concated_ints = self.tas.concated_ints.artifact
        iterations = random.randint(1, 10)
        for _ in range(iterations):
            next_concated_ints, = concat_ints(next_concated_ints,
                                              next_concated_ints,
                                              self.tas.int_seq2.artifact,
                                              4, 6)

        fp = os.path.join(self.tempdir, 'very-concated-ints.qza')
        next_concated_ints.save(fp)
        dag = ProvDAG(fp)
        # iterations + o.g. concated_ints + o.g. int_seq + o.g. int_seq2
        self.assertEqual(len(dag), iterations + 3)

    def test_get_outer_provenance_nodes(self):
        fp = os.path.join(self.tempdir, 'disconnected-provenances')
        os.mkdir(fp)
        self.tas.concated_ints.artifact.save(
            os.path.join(fp, 'concated-ints.qza'))
        unattached_int_seq = Artifact.import_data('IntSequence1', [8, 8])
        unattached_int_seq.save(os.path.join(fp, 'unattached-int-seq.qza'))
        dag = ProvDAG(fp)

        actual = dag.get_outer_provenance_nodes(self.tas.concated_ints.uuid)
        exp = {
            self.tas.concated_ints.uuid,
            self.tas.int_seq1.uuid,
            self.tas.int_seq2.uuid
        }
        self.assertEqual(actual, exp)

    def test_get_outer_provenance_nodes_pipeline(self):
        dag = self.tas.pipeline_viz.dag
        actual = dag.get_outer_provenance_nodes(self.tas.pipeline_viz.uuid)
        exp = {
            self.tas.pipeline_viz.uuid,
            self.tas.int_seq1.uuid,
            self.tas.mapping1.uuid
        }
        self.assertEqual(actual, exp)

    def test_relabel_nodes(self):
        dag = ProvDAG(self.tas.concated_ints.filepath)

        exp_nodes = ['Harry', 'Ron', 'Hermione']
        new_labels = {}
        for node, exp_node in zip(dag.nodes, exp_nodes):
            new_labels[node] = exp_node
        dag.relabel_nodes(new_labels)

        # have all nodes been relabeled as expected?
        for exp_node in exp_nodes:
            self.assertIn(exp_node, dag.nodes)

        # are the UUIDs stored in the ProvNode payloads updated correctly?
        for node in dag.nodes:
            self.assertEqual(node, dag.get_node_data(node)._uuid)

    def test_relabel_nodes_with_copy(self):
        dag = self.tas.concated_ints.dag
        exp_nodes = ['Harry', 'Ron', 'Hermione']
        new_labels = {}
        for node, exp_node in zip(dag.nodes, exp_nodes):
            new_labels[node] = exp_node
        copied_dag = dag.relabel_nodes(new_labels, copy=True)

        self.assertIsInstance(copied_dag, ProvDAG)

        # have all nodes been relabeled as expected?
        for exp_node in exp_nodes:
            self.assertIn(exp_node, copied_dag.nodes)

        # are the UUIDs stored in the ProvNode payloads updated correctly?
        for node in copied_dag.nodes:
            self.assertEqual(node, copied_dag.get_node_data(node)._uuid)

    def test_collapsed_view(self):
        view = self.tas.concated_ints.dag.collapsed_view
        self.assertIsInstance(view, DiGraph)
        self.assertEqual(len(view), 3)

        exp_nodes = [
            self.tas.concated_ints.uuid,
            self.tas.int_seq1.uuid,
            self.tas.int_seq2.uuid
        ]
        for exp_node in exp_nodes:
            self.assertIn(exp_node, view.nodes)

    def test_collapsed_view_pipeline(self):
        view = self.tas.pipeline_viz.dag.collapsed_view
        self.assertIsInstance(view, DiGraph)
        self.assertEqual(len(view), 3)

        exp_nodes = [
            self.tas.pipeline_viz.uuid,
            self.tas.int_seq1.uuid,
            self.tas.mapping1.uuid
        ]
        for exp_node in exp_nodes:
            self.assertIn(exp_node, view.nodes)

    def test_invalid_provenance(self):
        """
        Mangle an intact v5 Archive so that its checksums.md5 is invalid,
        and then build a ProvDAG with it to confirm the ProvDAG constructor
        handles broken checksums appropriately
        """
        uuid = self.tas.int_seq1.uuid
        with generate_archive_with_file_removed(
            self.tas.int_seq1.filepath,
            uuid,
            os.path.join('data', 'ints.txt')
        ) as altered_archive:
            new_fp = os.path.join(uuid, 'data', 'tamper.txt')
            overwrite_fp = os.path.join(uuid, 'provenance', 'citations.bib')
            with zipfile.ZipFile(altered_archive, 'a') as zf:
                zf.writestr(new_fp, 'added file')

                with zf.open(overwrite_fp, 'w') as fh:
                    fh.write(b'999\n')

            expected = (
                '(?s)'
                f'Checksums are invalid for Archive {uuid}.*'
                'Archive may be corrupt.*'
                'Files added.*tamper.txt.*'
                'Files removed.*ints.txt.*'
                'Files changed.*provenance.*citations.bib.*'
            )

            with self.assertWarnsRegex(UserWarning, expected):
                dag = ProvDAG(altered_archive)

            self.assertEqual(dag.provenance_is_valid, ValidationCode.INVALID)

            diff = dag.checksum_diff
            self.assertEqual(list(diff.removed.keys()), ['data/ints.txt'])
            self.assertEqual(list(diff.added.keys()), ['data/tamper.txt'])
            self.assertEqual(list(diff.changed.keys()),
                             ['provenance/citations.bib'])

    def test_missing_checksums_md5(self):
        uuid = self.tas.single_int.uuid
        with generate_archive_with_file_removed(
            self.tas.single_int.filepath,
            uuid,
            'checksums.md5'
        ) as altered_archive:
            expected = (f'no item.*{uuid}.*Archive may be corrupt')
            with self.assertWarnsRegex(UserWarning, expected):
                dag = ProvDAG(altered_archive)

            self.assertEqual(dag.provenance_is_valid, ValidationCode.INVALID)

            diff = dag.checksum_diff
            self.assertEqual(diff, None)

    def test_error_if_missing_node_files(self):
        path_prefix = os.path.join('provenance', 'artifacts')
        root_uuid = self.tas.concated_ints.uuid
        for removed_file in [
            'metadata.yaml',
            'citations.bib',
            'VERSION',
            'action/action.yaml'
        ]:
            for uuid in [self.tas.int_seq1.uuid, self.tas.int_seq2.uuid]:
                with generate_archive_with_file_removed(
                    self.tas.concated_ints.filepath,
                    root_uuid,
                    os.path.join(path_prefix, uuid, removed_file)
                ) as altered_archive:
                    if removed_file == 'action/action.yaml':
                        file = 'action.yaml'
                    else:
                        file = removed_file

                    expected = (
                        f'(?s)Malformed.*{file}.*{uuid}.*corrupt.*'
                    )
                    with self.assertRaisesRegex(ValueError, expected):
                        ProvDAG(altered_archive)

    def test_v0_archive(self):
        dag = self.tas.table_v0.dag
        uuid = self.tas.table_v0.uuid

        self.assertEqual(dag.provenance_is_valid,
                         ValidationCode.PREDATES_CHECKSUMS)
        self.assertEqual(dag.node_has_provenance(uuid), False)

    def test_v1_archive(self):
        dag = self.tas.concated_ints_v1.dag
        uuid = self.tas.concated_ints_v1.uuid

        self.assertEqual(dag.provenance_is_valid,
                         ValidationCode.PREDATES_CHECKSUMS)
        self.assertEqual(dag.node_has_provenance(uuid), True)

    def test_v4_archive(self):
        dag = self.tas.concated_ints_v4.dag
        uuid = self.tas.concated_ints_v4.uuid

        self.assertEqual(dag.provenance_is_valid,
                         ValidationCode.PREDATES_CHECKSUMS)
        self.assertEqual(dag.node_has_provenance(uuid), True)

        with zipfile.ZipFile(self.tas.concated_ints_v4.filepath) as zf:
            citations_path = os.path.join(uuid, 'provenance', 'citations.bib')
            self.assertIn(citations_path, zf.namelist())

    def test_v5_archive(self):
        dag = self.tas.concated_ints_v5.dag
        uuid = self.tas.concated_ints_v5.uuid

        self.assertEqual(dag.provenance_is_valid,
                         ValidationCode.VALID)
        self.assertEqual(dag.node_has_provenance(uuid), True)

    def test_artifact_passed_as_metadata_archive(self):
        dag = self.tas.mapping1.dag
        uuid = self.tas.mapping1.uuid
        self.assertEqual(dag.node_has_provenance(uuid), True)
        self.assertEqual(dag.get_node_data(uuid)._uuid, uuid)
        self.assertEqual(dag.get_node_data(uuid).type, 'Mapping')

    def test_artifact_with_collection_of_inputs(self):
        dag = self.tas.merged_mappings.dag
        uuid = self.tas.merged_mappings.uuid
        root_node = dag.get_node_data(uuid)
        self.assertEqual(root_node.type, 'Mapping')

        exp_parents = {self.tas.mapping1.uuid, self.tas.mapping2.uuid}
        self.assertEqual(dag.predecessors(uuid), exp_parents)

    def test_provdag_initialized_from_provdag(self):
        for dag in [self.tas.single_int.dag, self.tas.concated_ints.dag,
                    self.tas.merged_mappings.dag]:
            copied = ProvDAG(dag)
            self.assertEqual(dag, copied)
            self.assertIsNot(dag, copied)

    def test_union_zero_or_one_dags(self):
        with self.assertRaisesRegex(ValueError, "pass.*two ProvDAGs"):
            ProvDAG.union([])

        with self.assertRaisesRegex(ValueError, "pass.*two ProvDAGs"):
            ProvDAG.union([self.tas.single_int.dag])

    def test_union_identity(self):
        dag = self.tas.single_int.dag
        uuid = self.tas.single_int.uuid
        unioned_dag = ProvDAG.union([dag, dag])

        self.assertEqual(dag, unioned_dag)
        self.assertSetEqual({uuid}, unioned_dag._parsed_artifact_uuids)
        self.assertEqual(unioned_dag.provenance_is_valid, ValidationCode.VALID)
        self.assertRegex(
            repr(unioned_dag),
            f'ProvDAG representing these Artifacts.*{uuid}')

    def test_union_two(self):
        unioned_dag = ProvDAG.union([self.tas.single_int.dag,
                                     self.tas.int_seq2.dag])

        self.assertEqual({self.tas.single_int.uuid, self.tas.int_seq2.uuid},
                         unioned_dag._parsed_artifact_uuids)
        self.assertEqual(unioned_dag.provenance_is_valid, ValidationCode.VALID)

        rep = repr(unioned_dag)
        self.assertRegex(rep, 'ProvDAG representing these Artifacts {')
        self.assertRegex(rep, f'{self.tas.single_int.uuid}')
        self.assertRegex(rep, f'{self.tas.int_seq2.uuid}')

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 2)

    def test_union_many(self):
        unioned_dag = ProvDAG.union([
            self.tas.single_int.dag,
            self.tas.int_seq1.dag,
            self.tas.mapping1.dag
        ])
        self.assertEqual(
            {
                self.tas.single_int.uuid,
                self.tas.int_seq1.uuid,
                self.tas.mapping1.uuid
            },
            unioned_dag._parsed_artifact_uuids
        )
        self.assertEqual(unioned_dag.provenance_is_valid, ValidationCode.VALID)

        rep = repr(unioned_dag)
        self.assertRegex(rep, 'ProvDAG representing these Artifacts {')
        self.assertRegex(rep, f'{self.tas.single_int.uuid}')
        self.assertRegex(rep, f'{self.tas.int_seq1.uuid}')
        self.assertRegex(rep, f'{self.tas.mapping1.uuid}')

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 3)

    def test_union_self_missing_checksums_md5(self):
        unioned_dag = ProvDAG.union([self.dag_missing_md5,
                                    self.tas.single_int.dag])

        self.assertRegex(repr(unioned_dag),
                         'ProvDAG representing these Artifacts.*'
                         f'{self.tas.single_int.uuid}')

        # The ChecksumDiff==None from the tinkered dag gets ignored...
        self.assertEqual(unioned_dag.checksum_diff, ChecksumDiff({}, {}, {}))

        # ...but this should make clear that the provenance is bad
        # (or that the user opted out of validation)
        self.assertEqual(unioned_dag.provenance_is_valid,
                         ValidationCode.INVALID)

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1)

    def test_union_other_missing_checksums_md5(self):
        """
        Tests unions of v5 dags where the other ProvDAG is missing its
        checksums.md5 but the calling ProvDAG is not
        """
        unioned_dag = ProvDAG.union([self.tas.single_int.dag,
                                    self.dag_missing_md5])

        self.assertRegex(repr(unioned_dag),
                         'ProvDAG representing these Artifacts.*'
                         f'{self.tas.single_int.uuid}')

        self.assertEqual(unioned_dag.checksum_diff, ChecksumDiff({}, {}, {}))

        self.assertEqual(unioned_dag.provenance_is_valid,
                         ValidationCode.INVALID)

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1)

    def test_union_both_missing_checksums_md5(self):
        """
        Tests unions of v5 dags where both artifacts are missing their
        checksums.md5 files.
        """
        unioned_dag = ProvDAG.union([self.dag_missing_md5,
                                    self.dag_missing_md5])

        self.assertRegex(repr(unioned_dag),
                         'ProvDAG representing these Artifacts.*'
                         f'{self.tas.single_int.uuid}')

        # Both DAGs have NoneType checksum_diffs, so the ChecksumDiff==None
        self.assertEqual(unioned_dag.checksum_diff, None)

        self.assertEqual(unioned_dag.provenance_is_valid,
                         ValidationCode.INVALID)

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1)

    def test_union_v0_v1_archives(self):
        unioned_dag = ProvDAG.union([self.tas.table_v0.dag,
                                     self.tas.concated_ints_v1.dag])

        self.assertIn(f'{self.tas.table_v0.uuid}', repr(unioned_dag))
        self.assertIn(f'{self.tas.concated_ints_v1.uuid}', repr(unioned_dag))

        self.assertEqual(unioned_dag.provenance_is_valid,
                         ValidationCode.PREDATES_CHECKSUMS)

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 2)

        self.assertFalse(
            unioned_dag.node_has_provenance(self.tas.table_v0.uuid))
        self.assertTrue(
            unioned_dag.node_has_provenance(self.tas.concated_ints_v1.uuid))

    def test_union_v1_v5_archives(self):
        unioned_dag = ProvDAG.union([self.tas.concated_ints_v1.dag,
                                     self.tas.concated_ints_v5.dag])

        self.assertIn(f'{self.tas.concated_ints_v1.uuid}', repr(unioned_dag))
        self.assertIn(f'{self.tas.concated_ints_v5.uuid}', repr(unioned_dag))

        self.assertEqual(unioned_dag.provenance_is_valid,
                         ValidationCode.PREDATES_CHECKSUMS)

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 2)

        self.assertTrue(
            unioned_dag.node_has_provenance(self.tas.concated_ints_v1.uuid))
        self.assertTrue(
            unioned_dag.node_has_provenance(self.tas.concated_ints_v5.uuid))

    def test_union_v5_v6_archives(self):
        unioned_dag = ProvDAG.union([self.tas.concated_ints_v5.dag,
                                     self.tas.concated_ints_v6.dag])

        self.assertIn(f'{self.tas.concated_ints_v5.uuid}', repr(unioned_dag))
        self.assertIn(f'{self.tas.concated_ints_v6.uuid}', repr(unioned_dag))

        self.assertEqual(unioned_dag.provenance_is_valid,
                         ValidationCode.VALID)

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 2)

        self.assertTrue(
            unioned_dag.node_has_provenance(self.tas.concated_ints_v5.uuid))
        self.assertTrue(
            unioned_dag.node_has_provenance(self.tas.concated_ints_v6.uuid))

    def test_dag_is_superset(self):
        """
        Tests union of three dags, where one dag is a proper superset of the
        others. We expect three _parsed_artifact_uuids, one terminal uuid,
        and one weakly_connected_component.
        """
        unioned_dag = ProvDAG.union([
            self.tas.concated_ints.dag,
            self.tas.int_seq1.dag,
            self.tas.int_seq2.dag
        ])

        self.assertIn(self.tas.int_seq1.uuid,
                      unioned_dag._parsed_artifact_uuids)
        self.assertIn(self.tas.int_seq2.uuid,
                      unioned_dag._parsed_artifact_uuids)
        self.assertIn(self.tas.concated_ints.uuid,
                      unioned_dag._parsed_artifact_uuids)
        self.assertEqual(len(unioned_dag._parsed_artifact_uuids), 3)

        self.assertEqual(len(unioned_dag.terminal_uuids), 1)
        self.assertEqual(unioned_dag.terminal_uuids,
                         {self.tas.concated_ints.uuid})

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1)

        # == tests identity of objects in memory, so we need is_isomorphic
        self.assertTrue(nx.is_isomorphic(self.tas.concated_ints.dag.dag,
                                         unioned_dag.dag))
        self.assertEqual(self.tas.concated_ints.dag, unioned_dag)

    def test_three_artifacts_two_terminal_uuids(self):
        """
        Tests union of four dags, where two artifacts have shared parents but
        no direct relationship with eachother. We expect four
        _parsed_artifact_uuids, two terminal uuids, and one
        weakly_connected_component.
        """
        unioned_dag = ProvDAG.union([
            self.tas.int_seq1.dag,
            self.tas.int_seq2.dag,
            self.tas.concated_ints.dag,
            self.tas.other_concated_ints.dag
        ])

        self.assertIn(self.tas.int_seq1.uuid,
                      unioned_dag._parsed_artifact_uuids)
        self.assertIn(self.tas.int_seq2.uuid,
                      unioned_dag._parsed_artifact_uuids)
        self.assertIn(self.tas.concated_ints.uuid,
                      unioned_dag._parsed_artifact_uuids)
        self.assertIn(self.tas.other_concated_ints.uuid,
                      unioned_dag._parsed_artifact_uuids)
        self.assertEqual(len(unioned_dag._parsed_artifact_uuids), 4)

        self.assertEqual(len(unioned_dag.terminal_uuids), 2)
        self.assertEqual(
            unioned_dag.terminal_uuids,
            {self.tas.concated_ints.uuid, self.tas.other_concated_ints.uuid}
        )

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1)

    def test_one_analysis_two_artifacts(self):
        """
        In this set of test archives, both artifacts are derived from the same
        parents, so should produce one connected DAG even though we are
        missing the parent artifacts used to create them
        """
        unioned_dag = ProvDAG.union([
            self.tas.concated_ints.dag,
            self.tas.other_concated_ints.dag
        ])

        self.assertIn(self.tas.concated_ints.uuid,
                      unioned_dag._parsed_artifact_uuids)
        self.assertIn(self.tas.other_concated_ints.uuid,
                      unioned_dag._parsed_artifact_uuids)
        self.assertEqual(len(unioned_dag._parsed_artifact_uuids), 2)

        self.assertEqual(len(unioned_dag.terminal_uuids), 2)
        self.assertEqual(
            unioned_dag.terminal_uuids,
            {self.tas.concated_ints.uuid, self.tas.other_concated_ints.uuid}
        )

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1)

    def test_no_checksum_validation_on_intact_artifact(self):
        no_validation_dag = ProvDAG(self.tas.concated_ints.filepath,
                                    validate_checksums=False)

        self.assertEqual(len(no_validation_dag.terminal_uuids), 1)
        terminal_uuid, *_ = no_validation_dag.terminal_uuids
        self.assertEqual(terminal_uuid, self.tas.concated_ints.uuid)

        self.assertEqual(len(no_validation_dag), 3)
        self.assertEqual(no_validation_dag.provenance_is_valid,
                         ValidationCode.VALIDATION_OPTOUT)
        self.assertEqual(no_validation_dag.checksum_diff, None)

    def test_no_checksum_validation_missing_checksums_md5(self):
        with generate_archive_with_file_removed(
            self.tas.concated_ints.filepath,
            self.tas.concated_ints.uuid,
            'checksums.md5'
        ) as altered_archive:
            dag = ProvDAG(altered_archive, validate_checksums=False)

            self.assertEqual(
                dag.provenance_is_valid, ValidationCode.VALIDATION_OPTOUT)

            self.assertEqual(dag.checksum_diff, None)

    def test_no_checksum_validation_missing_node_files(self):
        path_prefix = os.path.join('provenance', 'artifacts')
        root_uuid = self.tas.concated_ints.uuid
        for removed_file in [
            'metadata.yaml',
            'citations.bib',
            'VERSION',
            'action/action.yaml'
        ]:
            for uuid in [self.tas.int_seq1.uuid, self.tas.int_seq2.uuid]:
                with generate_archive_with_file_removed(
                    self.tas.concated_ints.filepath,
                    root_uuid,
                    os.path.join(path_prefix, uuid, removed_file)
                ) as altered_archive:
                    if removed_file == 'action/action.yaml':
                        file = 'action.yaml'
                    else:
                        file = removed_file

                    expected = (
                        f'(?s)Malformed.*{file}.*{uuid}.*corrupt.*'
                    )
                    with self.assertRaisesRegex(ValueError, expected):
                        ProvDAG(altered_archive, validate_checksums=False)


class EmptyParserTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp(prefix='qiime2-test-parse-temp-')

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_get_parser(self):
        parser = EmptyParser.get_parser(None)
        self.assertIsInstance(parser, EmptyParser)

    def test_get_parser_input_data_not_none(self):
        fn = 'not-a-zip.txt'
        fp = os.path.join(self.tempdir, fn)
        with open(fp, 'w') as fh:
            fh.write('some text\n')
        with self.assertRaisesRegex(
                TypeError, f"EmptyParser.*{fn} is not None"):
            EmptyParser.get_parser(fp)

    def test_parse_a_nonetype(self):
        """
        tests that we can actually create empty ProvDAGs
        """
        parser = EmptyParser()
        parsed = parser.parse_prov(Config(), None)
        self.assertIsInstance(parsed, ParserResults)
        self.assertEqual(parsed.parsed_artifact_uuids, set())
        self.assertTrue(nx.is_isomorphic(parsed.prov_digraph, nx.DiGraph()))
        self.assertEqual(parsed.provenance_is_valid, ValidationCode.VALID)
        self.assertEqual(parsed.checksum_diff, None)


class ProvDAGParserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dags = dict()
        for archive_version in TEST_DATA:
            # supress warning from parsing provenance for a v0 provDag
            uuid = TEST_DATA['0']['uuid']
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore',  f'Art.*{uuid}.*prior')
                cls.dags[archive_version] = ProvDAG(
                    str(TEST_DATA[archive_version]['qzv_fp']))

    def setUp(self):
        self.tempdir = tempfile.mkdtemp(prefix='qiime2-test-parse-temp-')

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    # TODO: make artifacts of all versions, then test
    def test_get_parser(self):
        for version in TEST_DATA:
            parser = ProvDAGParser.get_parser(self.dags[version])
            self.assertIsInstance(parser, ProvDAGParser)

    def test_get_parser_input_data_not_a_provdag(self):
        fn = 'not_a_zip.txt'
        fp = os.path.join(DATA_DIR, fn)
        with self.assertRaisesRegex(
                TypeError, f"ProvDAGParser.*{fn} is not a ProvDAG"):
            ProvDAGParser.get_parser(fp)

    def test_parse_a_provdag(self):
        parser = ProvDAGParser()
        for dag in self.dags.values():
            parsed = parser.parse_prov(Config(), dag)
            self.assertIsInstance(parsed, ParserResults)
            self.assertEqual(parsed.parsed_artifact_uuids,
                             dag._parsed_artifact_uuids)
            # NOTE: networkx thinks about graph equality in terms of object
            # identity, so we must use nx.is_isomorphic to confirm "equality"
            self.assertTrue(nx.is_isomorphic(parsed.prov_digraph, dag.dag))
            self.assertEqual(parsed.provenance_is_valid,
                             dag.provenance_is_valid)
            self.assertEqual(parsed.checksum_diff, dag.checksum_diff)


class SelectParserTests(unittest.TestCase):
    def test_correct_parser_type(self):
        empty = select_parser(None)
        self.assertIsInstance(empty, EmptyParser)

        test_arch_fp = TEST_DATA['5']['qzv_fp']
        archive = select_parser(test_arch_fp)
        self.assertIsInstance(archive, ArchiveParser)

        dag = ProvDAG()
        pdag = select_parser(dag)
        self.assertIsInstance(pdag, ProvDAGParser)

        # check dir_fp as fp
        dir_fp = pathlib.Path(DATA_DIR) / 'parse_dir_test'
        dir_p = select_parser(dir_fp)
        self.assertIsInstance(dir_p, DirectoryParser)

        # check dir_fp as str
        dir_fp_str = str(dir_fp)
        dir_p = select_parser(dir_fp_str)
        self.assertIsInstance(dir_p, DirectoryParser)

    def test_correct_archive_parser_version(self):
        for arch_ver in TEST_DATA:
            qzv_fp = TEST_DATA[arch_ver]['qzv_fp']
            handler = select_parser(qzv_fp)
            self.assertIsInstance(handler, TEST_DATA[arch_ver]['parser'])


class ParseProvenanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = Config()

    def test_parse_with_artifact_parser(self):
        uuid = TEST_DATA['5']['uuid']
        qzv_fp = TEST_DATA['5']['qzv_fp']
        parser_results = parse_provenance(self.cfg, qzv_fp)
        self.assertIsInstance(parser_results, ParserResults)
        p_a_uuids = parser_results.parsed_artifact_uuids
        self.assertIsInstance(p_a_uuids, set)
        self.assertIsInstance(next(iter(p_a_uuids)), UUID)
        self.assertEqual(len(parser_results.prov_digraph), 15)
        self.assertIn(uuid, parser_results.prov_digraph)
        self.assertIsInstance(
            parser_results.prov_digraph.nodes[uuid]['node_data'],
            ProvNode)
        self.assertEqual(parser_results.provenance_is_valid,
                         TEST_DATA['5']['prov_is_valid'])
        self.assertEqual(parser_results.checksum_diff,
                         TEST_DATA['5']['checksum'])

    def test_parse_with_provdag_parser(self):
        uuid = TEST_DATA['5']['uuid']
        qzv_fp = TEST_DATA['5']['qzv_fp']
        starter = ProvDAG(qzv_fp)
        parser_results = parse_provenance(self.cfg, starter)
        self.assertIsInstance(parser_results, ParserResults)
        p_a_uuids = parser_results.parsed_artifact_uuids
        self.assertIsInstance(p_a_uuids, set)
        self.assertIsInstance(next(iter(p_a_uuids)), UUID)
        self.assertEqual(len(parser_results.prov_digraph), 15)
        self.assertIn(uuid, parser_results.prov_digraph)
        self.assertIsInstance(
            parser_results.prov_digraph.nodes[uuid]['node_data'],
            ProvNode)
        self.assertEqual(parser_results.provenance_is_valid,
                         TEST_DATA['5']['prov_is_valid'])
        self.assertEqual(parser_results.checksum_diff,
                         TEST_DATA['5']['checksum'])

    def test_parse_with_empty_parser(self):
        res = parse_provenance(self.cfg, None)
        self.assertIsInstance(res, ParserResults)
        self.assertEqual(res.parsed_artifact_uuids, set())
        self.assertTrue(
            nx.is_isomorphic(res.prov_digraph, nx.DiGraph()))
        self.assertEqual(res.provenance_is_valid, ValidationCode.VALID)
        self.assertEqual(res.checksum_diff, None)

    def test_parse_with_directory_parser(self):
        # Non-recursive
        dir_fp = pathlib.Path(DATA_DIR) / 'parse_dir_test'
        res = parse_provenance(self.cfg, dir_fp)
        self.assertEqual(self.cfg.recurse, False)
        self.assertIsInstance(res, ParserResults)
        v5_uu_id = 'ffb7cee3-2f1f-4988-90cc-efd5184ef003'
        self.assertEqual(res.parsed_artifact_uuids, {v5_uu_id})
        self.assertEqual(len(res.prov_digraph), 15)
        self.assertEqual(res.provenance_is_valid,
                         TEST_DATA['5']['prov_is_valid'])
        self.assertEqual(res.checksum_diff, TEST_DATA['5']['checksum'])

        # Recursive
        self.cfg.recurse = True
        self.assertEqual(self.cfg.recurse, True)
        res = parse_provenance(self.cfg, dir_fp)
        self.assertIsInstance(res, ParserResults)
        v5_unr_tree_id = '12e012d5-b01c-40b7-b825-a17f0478a02f'
        v5_tbl_id = '89af91c0-033d-4e30-8ac4-f29a3b407dc1'
        v5_uu_id = 'ffb7cee3-2f1f-4988-90cc-efd5184ef003'
        self.assertEqual(res.parsed_artifact_uuids,
                         {v5_tbl_id, v5_unr_tree_id, v5_uu_id})
        self.assertEqual(len(res.prov_digraph), 16)
        self.assertEqual(res.provenance_is_valid,
                         TEST_DATA['5']['prov_is_valid'])
        self.assertEqual(res.checksum_diff, TEST_DATA['5']['checksum'])

    def test_parse_with_directory_parser_bad_dir_path(self):
        dir_fp = pathlib.Path(DATA_DIR) / 'fake_dir'
        with self.assertRaisesRegex(UnparseableDataError, 'not a valid dir'):
            parse_provenance(self.cfg, dir_fp)

    def test_no_correct_parser_found_error(self):
        input_data = {'this': 'is not parseable'}
        with self.assertRaisesRegex(
            UnparseableDataError,
            f"(?s)Input data {input_data}.*not supported.*"
            "AttributeError.*ArchiveParser.*dict.*no attribute.*seek.*"
            "DirectoryParser.*expects a directory.*"
            "ProvDAGParser.*is not a ProvDAG.*"
            "EmptyParser.*is not None"
                ):
            select_parser(input_data)


class DirectoryParserTests(unittest.TestCase):
    def test_parse_empty_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError,
                                        f"No .qza or .qzv files.*{tmpdir}"):
                ProvDAG(tmpdir)

    def test_directory_parser_works_regardless_trailing_slash(self):
        dag = ProvDAG(DATA_DIR + '/parse_dir_test/')
        dag2 = ProvDAG(DATA_DIR + '/parse_dir_test')
        self.assertEqual(dag, dag2)

    def test_directory_parser_captures_all_parsed_artifact_uuids(self):
        """
        This test dir contains a feature table, an unrooted tree, and a viz.
        The feature table is a true subset of the viz, and is not re-parsed as
        implemented.

        The unrooted tree is a terminal node, but shares some ancestors with
        the others. It must be parsed, though its shared parents don't.

        The resulting dag should be aware of all three user-passed artifacts,
        regardless of whether we actually parse the v5_table
        """
        dag = ProvDAG(DATA_DIR + '/parse_dir_test/', recurse=True)
        v5_tbl_id = '89af91c0-033d-4e30-8ac4-f29a3b407dc1'
        v5_uu_id = 'ffb7cee3-2f1f-4988-90cc-efd5184ef003'
        v5_unr_tree_id = '12e012d5-b01c-40b7-b825-a17f0478a02f'
        self.assertEqual(dag._parsed_artifact_uuids,
                         {v5_tbl_id, v5_uu_id, v5_unr_tree_id})

    def test_archive_not_parsed(self):
        mixed_archive_fp = os.path.join(DATA_DIR, 'mixed_v0_v1_uu_emperor.qzv')
        v1_uuid = TEST_DATA['1']['uuid']
        v0_uuid = '9f6a0f3e-22e6-4c39-8733-4e672919bbc7'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'):
            dag = ProvDAG(mixed_archive_fp)
        # Never parsed
        self.assertTrue(archive_not_parsed(root_id='not even an id', dag=dag))
        # Only parsed as a !no-provenance parent
        self.assertTrue(archive_not_parsed(v0_uuid, dag))
        # Actually parsed
        self.assertFalse(archive_not_parsed(v1_uuid, dag))

    def test_directory_parser_idempotent_with_parse_and_union(self):
        # Non-recursive
        base_dir = os.path.join(DATA_DIR, 'parse_dir_test')
        inner_dir = os.path.join(base_dir, 'inner')
        dir_dag = ProvDAG(inner_dir)
        # parse files separately, then union
        tbl = ProvDAG(os.path.join(inner_dir, 'v5_table.qza'))
        tree = ProvDAG(os.path.join(inner_dir, 'v5_unrooted_tree.qza'))
        union_dag = ProvDAG.union([tbl, tree])
        self.assertEqual(dir_dag, union_dag)

        # Recursive
        base_dir = os.path.join(DATA_DIR, 'parse_dir_test')
        inner_dir = os.path.join(base_dir, 'inner')
        dir_dag = ProvDAG(base_dir, recurse=True)
        # parse files separately, then union
        tbl = ProvDAG(os.path.join(inner_dir, 'v5_table.qza'))
        tree = ProvDAG(os.path.join(inner_dir, 'v5_unrooted_tree.qza'))
        viz = ProvDAG(os.path.join(base_dir, 'v5_uu_emperor.qzv'))
        union_dag = ProvDAG.union([tbl, tree, viz])
        self.assertEqual(dir_dag, union_dag)

    def test_directory_parser_multiple_imports(self):
        base_dir = os.path.join(DATA_DIR, 'multiple_imports_test')
        inner_dir = os.path.join(base_dir, 'duplicated_inner')
        inner_dir_dag = ProvDAG(inner_dir)
        s1_id = '4f6794e7-0e34-46d9-9a48-3fbc7900430e'
        s2_id = 'b4fd43fb-91c3-45f6-9672-7cf8fd90bc0b'
        self.assertEqual(len(inner_dir_dag), 2)
        self.assertIn(s1_id, inner_dir_dag.dag)
        self.assertIn(s2_id, inner_dir_dag.dag)

        # Despite the two pairs of duplicate files with different names,
        # this DAG should be identical to the inner.
        dir_dag = ProvDAG(base_dir)
        self.assertEqual(len(inner_dir_dag), 2)
        self.assertIn(s1_id, inner_dir_dag.dag)
        self.assertIn(s2_id, inner_dir_dag.dag)

        self.assertEqual(dir_dag, inner_dir_dag)

    def test_verbose(self):
        buffer = io.StringIO()
        tbl = 'parse_dir_test/inner/v5_table.qza'
        viz = 'parse_dir_test/v5_uu_emperor.qzv'
        tree = 'parse_dir_test/inner/v5_unrooted_tree.qza'
        with redirect_stdout(buffer):
            dag = ProvDAG(os.path.join(DATA_DIR, 'parse_dir_test'),
                          verbose=True, recurse=True)
        self.assertEqual(dag.cfg.verbose, True)
        stdout_log = buffer.getvalue()
        self.assertRegex(stdout_log, f"parsing.*{tbl}")
        self.assertRegex(stdout_log, f"parsing.*{viz}")
        self.assertRegex(stdout_log, f"parsing.*{tree}")
