from contextlib import redirect_stdout
import copy
import io
import networkx as nx
import os
import pathlib
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


class ProvDAGTests(unittest.TestCase):
    def setUp(self):
        self.pm = PluginManager()
        self.dp = self.pm.plugins['dummy-plugin']
        self.tempdir = tempfile.mkdtemp(prefix='qiime2-test-parse-temp-')

        self.int_seq = Artifact.import_data('IntSequence1', [1, 2, 3])
        self.int_seq2 = Artifact.import_data('IntSequence2', [20, 30])
        concat_ints = self.dp.actions['concatenate_ints']
        self.concated_ints, = concat_ints(self.int_seq, self.int_seq,
                                          self.int_seq2, 9, 0)

    def tearDown(self):
        shutil.rmtree(self.tempdir)

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

    def test_number_of_nodes(self):
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        self.int_seq.save(fp)
        dag = ProvDAG(fp)
        self.assertEqual(len(dag.nodes), 1)

        fp = os.path.join(self.tempdir, 'concated_ints.qza')
        self.concated_ints.save(fp)
        dag = ProvDAG(fp)
        self.assertEqual(len(dag.nodes), 3)

    def test_number_of_terminal_nodes(self):
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        self.int_seq.save(fp)
        dag = ProvDAG(fp)
        self.assertEqual(len(dag.terminal_nodes), 1)

        fp = os.path.join(self.tempdir, 'concated_ints.qza')
        self.concated_ints.save(fp)
        dag = ProvDAG(fp)
        self.assertEqual(len(dag.terminal_nodes), 1)

        fp = os.path.join(self.tempdir, 'two-ints')
        os.mkdir(fp)
        self.int_seq.save(os.path.join(fp, 'int-seq1.qza'))
        self.int_seq2.save(os.path.join(fp, 'int-seq2.qza'))
        dag = ProvDAG(fp)
        self.assertEqual(len(dag.terminal_nodes), 2)

        fp = os.path.join(self.tempdir, 'int-and-concated-ints')
        os.mkdir(fp)
        self.int_seq.save(os.path.join(fp, 'int-seq1.qza'))
        self.concated_ints.save(os.path.join(fp, 'concated-ints.qza'))
        dag = ProvDAG(fp)
        self.assertEqual(len(dag.terminal_nodes), 1)

    def test_root_node_is_archive_root(self):
        fp = os.path.join(self.tempdir, 'concated-ints.qza')
        self.concated_ints.save(fp)
        dag = ProvDAG(fp)
        with zipfile.ZipFile(fp) as zf:
            all_filenames = zf.namelist()
            root_filenames = filter(is_root_provnode_data, all_filenames)
            root_filepaths = [pathlib.Path(fp) for fp in root_filenames]
            exp_node = ProvNode(Config(), zf, root_filepaths)
            act_terminal_node, *_ = dag.terminal_nodes
            self.assertEqual(exp_node, act_terminal_node)

    def test_number_of_actions(self):
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        self.int_seq.save(fp)
        dag = ProvDAG(fp)
        self.assertEqual(dag.dag.number_of_edges(), 0)

        fp = os.path.join(self.tempdir, 'concated-ints.qza')
        self.concated_ints.save(fp)
        dag = ProvDAG(fp)
        self.assertEqual(dag.dag.number_of_edges(), 2)

    def test_nonexistent_fp(self):
        fp = os.path.join(self.tempdir, 'does-not-exist.qza')
        err_msg = f'FileNotFoundError.*ArchiveParser.*{fp}'
        with self.assertRaisesRegex(UnparseableDataError, err_msg):
            ProvDAG(fp)

    def test_insufficient_permissions(self):
        fp = os.path.join(self.tempdir, 'wrong-permissions.qza')
        self.int_seq.save(fp)
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
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        self.int_seq.save(fp)
        dag = ProvDAG(fp)
        self.assertIsInstance(dag.dag, DiGraph)

        fp = os.path.join(self.tempdir, 'concated-ints.qza')
        self.concated_ints.save(fp)
        dag = ProvDAG(fp)
        self.assertIsInstance(dag.dag, DiGraph)

    def test_dag_attributes(self):
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        self.int_seq.save(fp)
        dag = ProvDAG(fp)

        terminal_node, *_ = dag.terminal_nodes
        self.assertIsInstance(terminal_node, ProvNode)

        self.assertEqual(dag.provenance_is_valid, ValidationCode.VALID)

        empty_checksum_diff = ChecksumDiff(added={}, removed={}, changed={})
        self.assertEqual(dag.checksum_diff, empty_checksum_diff)

    def test_node_action_names(self):
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        self.int_seq.save(fp)
        dag = ProvDAG(fp)
        int_seq_node, *_ = dag.terminal_nodes
        self.assertEqual(int_seq_node.action.action_name, 'import')

        fp = os.path.join(self.tempdir, 'concated-ints.qza')
        self.concated_ints.save(fp)
        dag = ProvDAG(fp)
        concated_ints_node, *_ = dag.terminal_nodes
        self.assertEqual(concated_ints_node.action.action_name,
                         'concatenate_ints')

    def test_has_correct_edges(self):
        fp = os.path.join(self.tempdir, 'concated-ints.qza')
        self.concated_ints.save(fp)
        dag = ProvDAG(fp)

        for node in dag.nodes:
            prov_node = dag.get_node_data(node)
            if prov_node.action.action_name == 'concatenate_ints':
                concated_ints_node = prov_node
            elif prov_node.action.action_name == 'import':
                if prov_node.type == 'IntSequence1':
                    int_seq_node = prov_node
                elif prov_node.type == 'IntSequence2':
                    int_seq2_node = prov_node

        edges = dag.dag.edges
        self.assertIn((int_seq_node._uuid, concated_ints_node._uuid), edges)
        self.assertIn((int_seq2_node._uuid, concated_ints_node._uuid), edges)
        self.assertNotIn((int_seq_node._uuid, int_seq2_node._uuid), edges)
        self.assertNotIn((concated_ints_node._uuid, int_seq_node._uuid), edges)

    def test_dag_repr(self):
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        self.int_seq.save(fp)
        dag = ProvDAG(fp)
        int_seq_node, *_ = dag.terminal_nodes
        uuid = int_seq_node._uuid

        exp_repr = f'ProvDAG.*Artifacts.*{uuid}'
        self.assertRegex(repr(dag), exp_repr)

    def test_node_repr(self):
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        self.int_seq.save(fp)
        dag = ProvDAG(fp)
        int_seq_node, *_ = dag.terminal_nodes
        uuid = int_seq_node._uuid
        type_ = int_seq_node.type
        format_ = int_seq_node.format

        exp_repr = f'(?s)UUID.*{uuid}.*Type:.*{type_}.*Data Format:.*{format_}'
        self.assertRegex(repr(int_seq_node), exp_repr)

    def test_dag_eq(self):
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        self.int_seq.save(fp)
        dag = ProvDAG(fp)
        self.assertEqual(dag, dag)
        self.assertEqual(ProvDAG(fp), ProvDAG(fp))

        int_seq_dag = dag
        fp = os.path.join(self.tempdir, 'int-seq2.qza')
        self.int_seq2.save(fp)
        int_seq2_dag = ProvDAG(fp)
        # because they are isomorphic
        self.assertEqual(int_seq_dag, int_seq2_dag)

    def test_dag_not_eq(self):
        fp = os.path.join(self.tempdir, 'int-seq.qza')
        self.int_seq.save(fp)
        int_seq_dag = ProvDAG(fp)

        fp = os.path.join(self.tempdir, 'concated-ints.qza')
        self.concated_ints.save(fp)
        concated_ints_dag = ProvDAG(fp)

        self.assertNotEqual(int_seq_dag, concated_ints_dag)

    def test_v5_captures_full_history(self):
        nodes = self.dags['5'].nodes
        self.assertEqual(len(nodes), 15)
        node_list = ['ffb7cee3-2f1f-4988-90cc-efd5184ef003',
                     '0af08fa8-48b7-4c6a-83c6-e0f766156343',
                     '3b7d36ff-37ab-4ac2-958b-6a547d442bcf',
                     '7ecf8954-e49a-4605-992e-99fcee397935',
                     '9cc3281a-fefb-408e-8cf0-10637a06d84a',
                     '025e723d-b367-4812-820a-ae8bf8b80af4',
                     '83a80bfd-8954-4571-8fc7-ac9e8435156e',
                     '89af91c0-033d-4e30-8ac4-f29a3b407dc1',
                     '99fa3670-aa1a-45f6-ba8e-803c976a1163',
                     '430a6575-86b3-4cf6-b72e-0f7fce3ed342',
                     'a35830e1-4535-47c6-aa23-be295a57ee1c',
                     'aea3994b-0888-41c1-8e8c-69f6615d07cf',
                     'bce3d09b-e296-4f2b-9af4-834db6412429',
                     'd32a5ea6-1ca1-4635-b522-2253568ae35b',
                     'f20cecd6-9f82-4bde-a013-eb327612dc4d',
                     ]
        self.assertEqual(len(nodes), 15)
        self.assertEqual(set(nodes), set(node_list))

        # Terminal/alias node
        root_parents = [
            {'table': '89af91c0-033d-4e30-8ac4-f29a3b407dc1'},
            {'phylogeny': 'bce3d09b-e296-4f2b-9af4-834db6412429'}]
        self.assertEqual(nodes[node_list[0]]['node_data']._parents,
                         root_parents)
        # non-alias node
        n1_parents = [{'table': '89af91c0-033d-4e30-8ac4-f29a3b407dc1'},
                      ]
        self.assertEqual(nodes[node_list[1]]['node_data']._parents,
                         n1_parents)
        # some other nodes
        n2_parents = [{'tree': 'd32a5ea6-1ca1-4635-b522-2253568ae35b'},
                      ]
        self.assertEqual(nodes[node_list[2]]['node_data']._parents,
                         n2_parents)
        n3_parents = [{'demultiplexed_seqs':
                       '99fa3670-aa1a-45f6-ba8e-803c976a1163'}]
        self.assertEqual(nodes[node_list[3]]['node_data']._parents,
                         n3_parents)
        # import node
        n10_parents = []
        self.assertEqual(nodes[node_list[10]]['node_data']._parents,
                         n10_parents)

    def test_v5_get_outer_provenance_nodes(self):
        exp = {'ffb7cee3-2f1f-4988-90cc-efd5184ef003',
               'bce3d09b-e296-4f2b-9af4-834db6412429',
               '89af91c0-033d-4e30-8ac4-f29a3b407dc1',
               '7ecf8954-e49a-4605-992e-99fcee397935',
               '99fa3670-aa1a-45f6-ba8e-803c976a1163',
               'a35830e1-4535-47c6-aa23-be295a57ee1c',
               }
        root_uuid = TEST_DATA['5']['uuid']
        actual = self.dags['5'].get_outer_provenance_nodes(root_uuid)
        self.assertEqual(actual, exp)

    def test_v5_relabel_nodes(self):
        # This function modifies labels in place by default,
        # so create a local ProvDAG to protect our test data
        dag = ProvDAG(str(TEST_DATA['5']['qzv_fp']))
        # Test new node names
        exp_nodes = ['ffb7cee3',
                     '0af08fa8',
                     '3b7d36ff',
                     '7ecf8954',
                     '9cc3281a',
                     '025e723d',
                     '83a80bfd',
                     '89af91c0',
                     '99fa3670',
                     '430a6575',
                     'a35830e1',
                     'aea3994b',
                     'bce3d09b',
                     'd32a5ea6',
                     'f20cecd6',
                     ]
        new_labels = {node: node[:8] for node in dag.nodes}
        dag.relabel_nodes(new_labels)
        # Have all nodes been relabeled as expected?
        for node in exp_nodes:
            self.assertIn(node, dag.nodes)
        # Are the UUIDs stored in the ProvNode payloads updated correctly?
        for node in dag.nodes:
            self.assertEqual(node, dag.get_node_data(node)._uuid)
        self.assertEqual(dag._terminal_uuids, None)

        # Confirm terminal_uuids state is consistent with the relabeled node
        # names
        self.assertEqual(len(dag.terminal_uuids), 1)
        # This is deterministic because there is one uuid in the set:
        terminal_uuid, *_ = dag.terminal_uuids
        self.assertEqual(terminal_uuid, exp_nodes[0])

    def test_v5_relabel_nodes_with_copy(self):
        exp_nodes = ['ffb7cee3',
                     '0af08fa8',
                     '3b7d36ff',
                     '7ecf8954',
                     '9cc3281a',
                     '025e723d',
                     '83a80bfd',
                     '89af91c0',
                     '99fa3670',
                     '430a6575',
                     'a35830e1',
                     'aea3994b',
                     'bce3d09b',
                     'd32a5ea6',
                     'f20cecd6',
                     ]
        new_labels = {node: node[:8] for node in self.dags['5'].nodes}
        new_dag = self.dags['5'].relabel_nodes(new_labels, copy=True)
        # Have all nodes been relabeled as expected?
        for node in exp_nodes:
            self.assertIn(node, new_dag.nodes)
        # Are the UUIDs stored in the ProvNode payloads updated correctly?
        for node in new_dag.nodes:
            self.assertEqual(node, new_dag.get_node_data(node)._uuid)
        self.assertEqual(set(exp_nodes), set(new_dag.nodes))
        self.assertEqual(new_dag._terminal_uuids, None)

        # Confirm terminal_uuids state is consistent with the relabeled node
        # names
        self.assertEqual(len(new_dag.terminal_uuids), 1)
        # This is deterministic because there is one uuid in the set:
        terminal_uuid, *_ = new_dag.terminal_uuids
        self.assertEqual(terminal_uuid, exp_nodes[0])

    def test_v5_collapsed_view(self):
        exp_nodes = {'ffb7cee3-2f1f-4988-90cc-efd5184ef003',
                     'bce3d09b-e296-4f2b-9af4-834db6412429',
                     '89af91c0-033d-4e30-8ac4-f29a3b407dc1',
                     '7ecf8954-e49a-4605-992e-99fcee397935',
                     '99fa3670-aa1a-45f6-ba8e-803c976a1163',
                     'a35830e1-4535-47c6-aa23-be295a57ee1c',
                     }
        view = self.dags['5'].collapsed_view
        self.assertIsInstance(view, DiGraph)
        self.assertEqual(len(view), 6)
        for node in exp_nodes:
            self.assertIn(node, view.nodes)

    def test_invalid_provenance(self):
        """
        Mangle an intact v5 Archive so that its checksums.md5 is invalid,
        and then build a ProvDAG with it to confirm the ProvDAG constructor
        handles broken checksums appropriately

        Specifically:
        - remove the root `<uuid>/metadata.yaml`
        - add a new file called '<uuid>/tamper.txt`
        - overwrite `<uuid>/data/index.html` with '999\n'

        Modified from test_checksum_validator.test_checksums_mismatch
        """
        original_archive = TEST_DATA['5']['qzv_fp']
        drop_file = pathlib.Path('data') / 'emperor.html'
        root_uuid = TEST_DATA['5']['uuid']
        fp_pfx = pathlib.Path(root_uuid)
        with generate_archive_with_file_removed(
            qzv_fp=original_archive,
            root_uuid=root_uuid,
                file_to_drop=drop_file) as chopped_archive:

            # We'll also add a new file
            with zipfile.ZipFile(chopped_archive, 'a') as zf:
                new_fn = str(fp_pfx / 'tamper.txt')
                zf.writestr(new_fn, 'extra file')
                # and overwrite an existing file with junk
                extant_fn = str(fp_pfx / 'data' / 'index.html')

                # we expect a warning that we're overwriting the filename
                # this CM stops the warning from propagating up to stderr/out
                with self.assertWarnsRegex(UserWarning, 'Duplicate name'):
                    with zf.open(extant_fn, 'w') as myfile:
                        myfile.write(b'999\n')

            # Is our bad-checksums warning message correct?
            uuid = TEST_DATA['5']['uuid']
            expected = ('(?s)'
                        f'Checksums are invalid for Archive {uuid}.*\n'
                        'Archive may be corrupt.*\n'
                        'Files added.*tamper.*296583.*\n'
                        'Files removed.*emperor.*c42b3.*\n'
                        'Files changed.*data.*index.*065031.*f47bc3.*'
                        )
            with self.assertWarnsRegex(UserWarning, expected):
                a_dag = ProvDAG(chopped_archive)

            # Have we set provenance_is_valid correctly?
            self.assertEqual(a_dag.provenance_is_valid,
                             ValidationCode.INVALID)

            # Is the diff correct?
            diff = a_dag.checksum_diff
            self.assertEqual(list(diff.removed.keys()),
                             ['data/emperor.html'])
            self.assertEqual(
                diff.added,
                {'tamper.txt': '296583001b00d2b811b5871b19e0ad28'})
            self.assertEqual(
                diff.changed,
                {'data/index.html': ('065031e17943cd0780f197874c4f011e',
                                     'f47bc36040d5c7db08e4b3a457dcfbb2')
                 })

    def test_v5_archive_has_invalid_checksums(self):
        """
        Remove a file from an intact v5 Archive so that its checksums.md5 is
        invalid, and then build a ProvDAG with it to confirm the ProvDAG
        constructor handles broken checksums appropriately
        """
        drop_file = pathlib.Path('data') / 'index.html'
        with generate_archive_with_file_removed(
            qzv_fp=TEST_DATA['5']['qzv_fp'],
            root_uuid=TEST_DATA['5']['uuid'],
                file_to_drop=drop_file) as chopped_archive:

            # Is our bad-checksums warning message correct?
            uuid = TEST_DATA['5']['uuid']
            expected = (f'(?s)Checksums are invalid for Archive {uuid}.*')
            with self.assertWarnsRegex(UserWarning, expected):
                a_dag = ProvDAG(chopped_archive)

            # Have we set provenance_is_valid correctly?
            self.assertEqual(a_dag.provenance_is_valid,
                             ValidationCode.INVALID)

            # Is the diff correct?
            diff = a_dag.checksum_diff
            self.assertEqual(list(diff.removed.keys()),
                             ['data/index.html'])
            self.assertEqual(diff.added, {})
            self.assertEqual(diff.changed, {})

    def test_v5_with_missing_checksums_md5(self):
        drop_file = pathlib.Path('checksums.md5')
        with generate_archive_with_file_removed(
            qzv_fp=TEST_DATA['5']['qzv_fp'],
            root_uuid=TEST_DATA['5']['uuid'],
                file_to_drop=drop_file) as chopped_archive:

            # Is our bad-checksums warning message correct?
            uuid = TEST_DATA['5']['uuid']
            expected = (f'no item.*{uuid}.*Archive may be corrupt')
            with self.assertWarnsRegex(UserWarning, expected):
                a_dag = ProvDAG(chopped_archive)

            # Have we set provenance_is_valid correctly?
            self.assertEqual(a_dag.provenance_is_valid,
                             ValidationCode.INVALID)

            # Is the diff correct?
            diff = a_dag.checksum_diff
            self.assertEqual(diff, None)

    def test_ProvDAG_error_if_missing_node_files(self):
        pfx = 'provenance/artifacts/'
        for archive_version in TEST_DATA:
            # V0 doesn't have root nodes
            if archive_version == '0':
                continue
            root_uuid = TEST_DATA[archive_version]['uuid']
            node_uuid = TEST_DATA[archive_version]['nonroot_node']
            parser = TEST_DATA[archive_version]['parser']
            fnames = parser.expected_files_in_all_nodes
            for name in fnames:
                drop_file = pathlib.Path(pfx) / node_uuid / name
                with generate_archive_with_file_removed(
                    qzv_fp=TEST_DATA[archive_version]['qzv_fp'],
                    root_uuid=root_uuid,
                        file_to_drop=drop_file) as chopped_archive:

                    # Fudging this to match what the user sees - 'action.yaml'
                    if name == 'action/action.yaml':
                        name = 'action.yaml'
                    fn = 'mangled.qzv'
                    expected = (
                        f"(?s)Malformed.*{name}.*{node_uuid}.*"
                        f"{fn}.*corrupt"
                    )
                    with self.assertRaisesRegex(ValueError, expected):
                        # Only v5 warns on this, so an assert would be clunky
                        with warnings.catch_warnings():
                            warnings.filterwarnings(
                                'ignore',
                                f'Checksums.*invalid.*{root_uuid}',
                                UserWarning)
                            ProvDAG(chopped_archive)

    def test_mixed_v0_v1_archive(self):
        mixed_archive_fp = os.path.join(DATA_DIR, 'mixed_v0_v1_uu_emperor.qzv')
        v1_uuid = TEST_DATA['1']['uuid']
        v0_uuid = '9f6a0f3e-22e6-4c39-8733-4e672919bbc7'

        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'):
            dag = ProvDAG(mixed_archive_fp)
            self.assertEqual(dag.node_has_provenance(v1_uuid), True)
            self.assertEqual(dag.get_node_data(v1_uuid)._uuid, v1_uuid)

            self.assertEqual(dag.node_has_provenance(v0_uuid), False)
            self.assertEqual(dag.get_node_data(v0_uuid), None)

    def test_artifact_passed_as_metadata_archive(self):
        """
        Tests:
        - smoke
        - does the parser find the captured provenance?
        - is the UUID parsed correctly?
        - is the node's type correct?
        """
        a_as_md_fp = os.path.join(DATA_DIR, 'artifact_as_md_v5.qzv')
        a_as_md_uuid = 'd1d36ada-29a5-436e-9136-304a8b25ff10'

        dag = ProvDAG(a_as_md_fp)
        self.assertEqual(dag.node_has_provenance(a_as_md_uuid), True)
        self.assertEqual(dag.get_node_data(a_as_md_uuid)._uuid, a_as_md_uuid)
        self.assertEqual(dag.get_node_data(a_as_md_uuid).type,
                         'FeatureData[Taxonomy]')

    def test_artifact_with_collection_of_inputs(self):
        fp = os.path.join(DATA_DIR, 'merged_tbls.qza')
        root_uuid = '2a045e27-7f3a-4d83-b358-8d39373708cb'
        dag = ProvDAG(fp)
        root_node = dag.get_node_data(root_uuid)
        self.assertEqual(root_node.type, 'FeatureTable[Frequency]')
        exp_parents = {
            '84898e39-f6e0-44bb-8fa1-6df2f330af68',
            '0be6c7be-ad84-4417-9f1c-cade0a8a9b58'
        }
        self.assertEqual(dag.predecessors(root_uuid), exp_parents)

    def test_provdag_initialized_from_a_provdag(self):
        for dag in self.dags.values():
            copied = ProvDAG(dag)
            self.assertEqual(dag, copied)
            self.assertIsNot(dag, copied)

    def test_predecessors(self):
        exp = {'bce3d09b-e296-4f2b-9af4-834db6412429',
               '89af91c0-033d-4e30-8ac4-f29a3b407dc1'}
        act = self.dags['5'].predecessors(TEST_DATA['5']['uuid'])
        self.assertSetEqual(exp, act)

    def test_predecessors_not_collapsed(self):
        inner_node = '83a80bfd-8954-4571-8fc7-ac9e8435156e'
        exp = {'9cc3281a-fefb-408e-8cf0-10637a06d84a'}
        this_dag = self.dags['5']
        act = this_dag.predecessors(inner_node, this_dag)
        self.assertSetEqual(exp, act)


class ProvDAGUnionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Because union is copy-only, we can create our test data once here,
        and union it to our hearts' content.
        """
        cls.v3_dag = ProvDAG(str(TEST_DATA['3']['qzv_fp']))
        cls.v4_dag = ProvDAG(str(TEST_DATA['4']['qzv_fp']))
        cls.v5_qzv = ProvDAG(str(TEST_DATA['5']['qzv_fp']))
        cls.v5_table = ProvDAG(os.path.join(DATA_DIR, 'v5_table.qza'))

        cls.v3_uuid = TEST_DATA['3']['uuid']
        cls.v4_uuid = TEST_DATA['4']['uuid']
        cls.qzv_uuid = TEST_DATA['5']['uuid']
        cls.table_uuid = '89af91c0-033d-4e30-8ac4-f29a3b407dc1'

        drop_file = pathlib.Path('checksums.md5')
        with generate_archive_with_file_removed(
            qzv_fp=TEST_DATA['5']['qzv_fp'],
            root_uuid=TEST_DATA['5']['uuid'],
                file_to_drop=drop_file) as chopped_archive:
            # we can't assert in a classmethod, so capture all warnings and
            # assert in test_setup_warnings
            with warnings.catch_warnings(record=True) as cls.w:
                cls.bad_dag = ProvDAG(chopped_archive)

    def test_setup_warnings(self):
        self.assertEqual(len(self.w), 1)
        with self.assertWarnsRegex(
            UserWarning,
                'There is no item named.*checksums.*in the archive'):
            warnings.warn(next(iter(self.w)))

    def test_union_zero_or_one_dags(self):
        """
        Tests union of zero or one ProvDAGs.
        """
        with self.assertRaisesRegex(ValueError, "pass.*two ProvDAGs"):
            ProvDAG.union([])

        with self.assertRaisesRegex(ValueError, "pass.*two ProvDAGs"):
            ProvDAG.union([self.v5_qzv])

    def test_union_identity(self):
        """
        Tests union of a ProvDAG with itself.
        """
        unioned_dag = ProvDAG.union([self.v5_qzv, self.v5_qzv])

        self.assertEqual(self.v5_qzv, unioned_dag)
        self.assertSetEqual({self.qzv_uuid},
                            unioned_dag._parsed_artifact_uuids)
        self.assertEqual(unioned_dag.provenance_is_valid, ValidationCode.VALID)
        self.assertRegex(
            repr(unioned_dag),
            f'ProvDAG representing these Artifacts.*{self.qzv_uuid}')

    def test_union_two(self):
        """
        Tests union of dag with another dag.
        Also checks that provenance_is_valid retains the lesser of the
        ValidationCodes from the v4- and v5+ dags.
        """
        unioned_dag = ProvDAG.union([self.v4_dag, self.v5_qzv])

        self.assertSetEqual({self.v4_uuid, self.qzv_uuid},
                            unioned_dag._parsed_artifact_uuids)
        self.assertEqual(unioned_dag.provenance_is_valid,
                         ValidationCode.PREDATES_CHECKSUMS)
        rep = repr(unioned_dag)
        self.assertRegex(rep, 'ProvDAG representing these Artifacts {')
        self.assertRegex(rep, f'{self.v4_uuid}')
        self.assertRegex(rep, f'{self.qzv_uuid}')

        # There should be two disconnected trees
        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 2)

    def test_union_many(self):
        """
        Tests union of dag with multiple other dags.
        Also checks that we have the correct number of disconnected trees
        (these three dags come from unrelated analyses, so should be disjoint)
        """
        unioned_dag = ProvDAG.union([self.v5_qzv, self.v4_dag, self.v3_dag])

        self.assertSetEqual({self.v3_uuid, self.v4_uuid, self.qzv_uuid},
                            unioned_dag._parsed_artifact_uuids)
        self.assertEqual(unioned_dag.provenance_is_valid,
                         ValidationCode.PREDATES_CHECKSUMS)
        rep = repr(unioned_dag)
        self.assertRegex(rep, 'ProvDAG representing these Artifacts {')
        self.assertRegex(rep, f'{self.v3_uuid}')
        self.assertRegex(rep, f'{self.v4_uuid}')
        self.assertRegex(rep, f'{self.qzv_uuid}')

        # There should be three disconnected trees
        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 3)

    def test_union_self_missing_checksums_md5(self):
        """
        Tests unions of v5 dags where the calling ProvDAG is missing its
        checksums.md5 but the other is not
        """
        unioned_dag = ProvDAG.union([self.bad_dag, self.v5_qzv])

        self.assertRegex(repr(unioned_dag),
                         'ProvDAG representing these Artifacts.*'
                         f'{self.qzv_uuid}')

        # The ChecksumDiff==None from the tinkered dag gets ignored...
        self.assertEqual(unioned_dag.checksum_diff, ChecksumDiff({}, {}, {}))

        # ...but this should make clear that the provenance is bad
        # (or that the user opted out of validation).
        self.assertEqual(unioned_dag.provenance_is_valid,
                         ValidationCode.INVALID)

        # There should be one fully-connected tree
        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1)

    def test_union_other_missing_checksums_md5(self):
        """
        Tests unions of v5 dags where the other ProvDAG is missing its
        checksums.md5 but the calling ProvDAG is not
        """
        unioned_dag = ProvDAG.union([self.v5_qzv, self.bad_dag])

        self.assertRegex(repr(unioned_dag),
                         'ProvDAG representing these Artifacts.*'
                         f'{self.qzv_uuid}')

        # The ChecksumDiff==None from the tinkered dag gets ignored...
        self.assertEqual(unioned_dag.checksum_diff, ChecksumDiff({}, {}, {}))

        # ...but this should make clear that the provenance is bad
        # (or that the user opted out of validation).
        self.assertEqual(unioned_dag.provenance_is_valid,
                         ValidationCode.INVALID)

        # There should be one fully-connected tree
        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1)

    def test_union_both_missing_checksums_md5(self):
        """
        Tests unions of v5 dags where both artifacts are missing their
        checksums.md5 files.
        """
        unioned_dag = ProvDAG.union([self.bad_dag, self.bad_dag])

        self.assertRegex(repr(unioned_dag),
                         'ProvDAG representing these Artifacts.*'
                         f'{self.qzv_uuid}')

        # Both DAGs have NoneType checksum_diffs, so the ChecksumDiff==None
        self.assertEqual(unioned_dag.checksum_diff, None)
        self.assertEqual(unioned_dag.provenance_is_valid,
                         ValidationCode.INVALID)

        # There should be one fully-connected tree
        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1)

    def test_one_dag_is_true_superset(self):
        """
        Tests union of three dags, where one dag is a true superset of the
        others. We expect three _parsed_artifact_uuids, one terminal uuid,
        and one weakly_connected_component.
        """
        v5_tree = ProvDAG(os.path.join(DATA_DIR, 'v5_rooted_tree.qza'))
        tree_uuid = 'bce3d09b-e296-4f2b-9af4-834db6412429'
        unmodified_dag = copy.copy(self.v5_qzv.dag)
        unioned_dag = ProvDAG.union([self.v5_qzv, self.v5_table, v5_tree])

        self.assertIn(self.qzv_uuid, unioned_dag._parsed_artifact_uuids)
        self.assertIn(self.table_uuid, unioned_dag._parsed_artifact_uuids)
        self.assertIn(tree_uuid, unioned_dag._parsed_artifact_uuids)
        self.assertEqual(len(unioned_dag._parsed_artifact_uuids), 3)

        self.assertEqual(len(unioned_dag.terminal_uuids), 1)
        self.assertEqual(unioned_dag.terminal_uuids, {self.qzv_uuid})

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1)

        # G == H tests identity of objects in memory, so we need
        # is_isomorphic
        self.assertTrue(nx.is_isomorphic(unmodified_dag, unioned_dag.dag))
        self.assertEqual(self.v5_qzv, unioned_dag)

    def test_three_artifacts_two_terminal_uuids(self):
        """
        Tests union of three dags, where the v5_qzv is downstream of the table,
        but not downstream of the unrooted tree. We expect three
        _parsed_artifact_uuids, two terminal uuids, and one
        weakly_connected_component.
        """
        v5_unr_tree = ProvDAG(os.path.join(DATA_DIR, 'v5_unrooted_tree.qza'))
        tree_uuid = '12e012d5-b01c-40b7-b825-a17f0478a02f'

        unioned_dag = ProvDAG.union([self.v5_qzv, self.v5_table, v5_unr_tree])

        self.assertIn(self.qzv_uuid, unioned_dag._parsed_artifact_uuids)
        self.assertIn(self.table_uuid, unioned_dag._parsed_artifact_uuids)
        self.assertIn(tree_uuid, unioned_dag._parsed_artifact_uuids)
        self.assertEqual(len(unioned_dag._parsed_artifact_uuids), 3)

        self.assertEqual(len(unioned_dag.terminal_uuids), 2)
        self.assertEqual(unioned_dag.terminal_uuids,
                         set([self.qzv_uuid, tree_uuid]))

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1)

    def test_graphs_same_analysis_missing_artifacts(self):
        """
        In this set of test archives, both .qzvs are derived from the same
        feature table, so should produce one connected DAG even though we are
        missing the rarefied_table.qza used to create the rarefied_table.qzv
        """
        rar_qzv = ProvDAG(os.path.join(DATA_DIR, 'v5_rarefied_table.qzv'))
        rar_uuid = '79a0d2ea-ea01-40c0-a4a4-0beab7c1f244'

        unioned_dag = ProvDAG.union([self.v5_qzv, self.v5_table, rar_qzv])

        self.assertIn(self.qzv_uuid, unioned_dag._parsed_artifact_uuids)
        self.assertIn(self.table_uuid, unioned_dag._parsed_artifact_uuids)
        self.assertIn(rar_uuid, unioned_dag._parsed_artifact_uuids)
        self.assertEqual(len(unioned_dag._parsed_artifact_uuids), 3)

        self.assertEqual(len(unioned_dag.terminal_uuids), 2)
        self.assertEqual(unioned_dag.terminal_uuids, {self.qzv_uuid, rar_uuid})

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1)


class ProvDAGTestsNoChecksumValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.no_checksum_dags = dict()
        for archive_version in TEST_DATA:
            # supress warning from parsing provenance for a v0 provDag
            uuid = TEST_DATA['0']['uuid']
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore',  f'Art.*{uuid}.*prior')
                cls.no_checksum_dags[archive_version] = ProvDAG(
                    str(TEST_DATA[archive_version]['qzv_fp']),
                    validate_checksums=False)

    # This should only trigger if something fails in setup or above
    # e.g. if a ProvDag fails to initialize
    def test_smoke(self):
        self.assertTrue(True)

    def test_no_checksum_validation_on_intact_artifact(self):
        dags = self.no_checksum_dags
        for vz in dags:
            self.assertEqual(len(dags[vz].terminal_uuids), 1)
            # This is deterministic because there is one uuid in the set:
            terminal_uuid, *_ = dags[vz].terminal_uuids
            self.assertEqual(terminal_uuid,
                             TEST_DATA[vz]['uuid'])
            # Node count acts as a proxy test for completeness here
            self.assertEqual(len(dags[vz]),
                             TEST_DATA[vz]['n_res'])
            self.assertEqual(dags[vz].provenance_is_valid,
                             ValidationCode.VALIDATION_OPTOUT)
            self.assertEqual(dags[vz].checksum_diff, None)

    def test_no_checksum_missing_checksums_md5(self):
        drop_file = pathlib.Path('checksums.md5')
        with generate_archive_with_file_removed(
            qzv_fp=TEST_DATA['5']['qzv_fp'],
            root_uuid=TEST_DATA['5']['uuid'],
                file_to_drop=drop_file) as chopped_archive:

            a_dag = ProvDAG(chopped_archive, validate_checksums=False)

            # Have we set provenance_is_valid correctly?
            self.assertEqual(
                a_dag.provenance_is_valid, ValidationCode.VALIDATION_OPTOUT)

            # Is the diff correct?
            diff = a_dag.checksum_diff
            self.assertEqual(diff, None)

    def test_no_checksum_validation_missing_node_files(self):
        pfx = 'provenance/artifacts/'
        for archive_version in TEST_DATA:
            # V0 doesn't have root nodes
            if archive_version == '0':
                continue
            root_uuid = TEST_DATA[archive_version]['uuid']
            node_uuid = TEST_DATA[archive_version]['nonroot_node']
            parser = TEST_DATA[archive_version]['parser']
            fnames = parser.expected_files_in_all_nodes
            for name in fnames:
                drop_file = pathlib.Path(pfx) / node_uuid / name
                with generate_archive_with_file_removed(
                    qzv_fp=TEST_DATA[archive_version]['qzv_fp'],
                    root_uuid=root_uuid,
                        file_to_drop=drop_file) as chopped_archive:

                    # Fudging this to match what the user sees - 'action.yaml'
                    if name == 'action/action.yaml':
                        name = 'action.yaml'
                    fn = 'mangled.qzv'
                    expected = (
                        f"(?s)Malformed.*{name}.*{node_uuid}.*"
                        f"{fn}.*corrupt"
                    )
                    with self.assertRaisesRegex(ValueError, expected):
                        ProvDAG(chopped_archive, validate_checksums=False)


class EmptyParserTests(unittest.TestCase):
    def test_get_parser(self):
        parser = EmptyParser.get_parser(None)
        self.assertIsInstance(parser, EmptyParser)

    def test_get_parser_input_data_not_none(self):
        fn = 'not_a_zip.txt'
        fp = os.path.join(DATA_DIR, fn)
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
