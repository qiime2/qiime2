import io
import os
import pathlib
import random
import shutil
import tempfile
import unittest
import zipfile
from contextlib import redirect_stdout

import networkx as nx
from networkx import DiGraph

from .._checksum_validator import ValidationCode
from ..parse import (
    ProvDAG, DirectoryParser, EmptyParser, ProvDAGParser, select_parser,
    parse_provenance, UnparseableDataError
)
from ..archive_parser import (
    ParserV0, ParserV1, ParserV2, ParserV3, ParserV4, ParserV5, ParserV6,
    Config, ProvNode, ParserResults, ArchiveParser,
)

from .testing_utilities import (
    is_root_provnode_data, generate_archive_with_file_removed, TestArtifacts
)

from qiime2 import Artifact
from qiime2.core.archive.archiver import ChecksumDiff


class ProvDAGTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def save_artifact_to_dir(self, artifact, directory):
        dir_fp = os.path.join(self.tempdir, directory)
        try:
            os.mkdir(dir_fp)
        except FileExistsError:
            pass

        fp = os.path.join(dir_fp, f'{artifact.name}.qza')
        artifact.artifact.save(fp)
        return dir_fp

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
        err_msg = 'FileNotFoundError'
        with self.assertRaisesRegex(UnparseableDataError, err_msg):
            ProvDAG(fp)

    def test_insufficient_permissions(self):
        fp = os.path.join(self.tempdir, 'int-seq-1-permissions-copy.qza')
        self.tas.int_seq1.artifact.save(fp)

        os.chmod(fp, 0o000)
        err_msg = 'PermissionError.*Permission denied'
        with self.assertRaisesRegex(UnparseableDataError, err_msg):
            ProvDAG(fp)

        os.chmod(fp, 0o123)
        err_msg = 'PermissionError.*Permission denied'
        with self.assertRaisesRegex(UnparseableDataError, err_msg):
            ProvDAG(fp)

    def test_not_a_zip_file(self):
        fp = os.path.join(self.tempdir, 'not-a-zip.txt')
        with open(fp, 'w') as fh:
            fh.write("This is just a text file.")

        err_msg = 'zipfile.BadZipFile.*File is not a zip file'
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
        concat_ints = self.tas.dp.actions['concatenate_ints']

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
        '''
        Mangle an intact v5 Archive so that its checksums.md5 is invalid,
        and then build a ProvDAG with it to confirm the ProvDAG constructor
        handles broken checksums appropriately
        '''
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
            expected = (
                'The checksums.md5 file is missing from the archive.*'
                'Archive may be corrupt'
            )
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

        self.assertEqual(
            dag.provenance_is_valid, ValidationCode.PREDATES_CHECKSUMS
        )
        self.assertEqual(dag.node_has_provenance(uuid), False)

    def test_v1_archive(self):
        dag = self.tas.concated_ints_v1.dag
        uuid = self.tas.concated_ints_v1.uuid

        self.assertEqual(
            dag.provenance_is_valid, ValidationCode.PREDATES_CHECKSUMS
        )
        self.assertEqual(dag.node_has_provenance(uuid), False)

    def test_v2_archive(self):
        dag = self.tas.concated_ints_v2.dag
        uuid = self.tas.concated_ints_v2.uuid

        self.assertEqual(
            dag.provenance_is_valid, ValidationCode.PREDATES_CHECKSUMS
        )
        self.assertEqual(dag.node_has_provenance(uuid), True)

    def test_v4_archive(self):
        dag = self.tas.concated_ints_v4.dag
        uuid = self.tas.concated_ints_v4.uuid

        self.assertEqual(
            dag.provenance_is_valid, ValidationCode.PREDATES_CHECKSUMS
        )
        self.assertEqual(dag.node_has_provenance(uuid), True)

        with zipfile.ZipFile(self.tas.concated_ints_v4.filepath) as zf:
            citations_path = os.path.join(uuid, 'provenance', 'citations.bib')
            self.assertIn(citations_path, zf.namelist())

    def test_v5_archive(self):
        dag = self.tas.concated_ints_v5.dag
        uuid = self.tas.concated_ints_v5.uuid

        self.assertEqual(dag.provenance_is_valid, ValidationCode.VALID)
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
            f'ProvDAG representing the provenance.*Artifacts.*{uuid}'
        )

    def test_union_two(self):
        unioned_dag = ProvDAG.union(
            [self.tas.single_int.dag, self.tas.int_seq2.dag])

        self.assertEqual(
            {self.tas.single_int.uuid, self.tas.int_seq2.uuid},
            unioned_dag._parsed_artifact_uuids
        )
        self.assertEqual(unioned_dag.provenance_is_valid, ValidationCode.VALID)

        rep = repr(unioned_dag)
        self.assertRegex(
            rep, 'ProvDAG representing the provenance.*Artifacts.'
        )
        self.assertRegex(rep, f'{self.tas.single_int.uuid}')
        self.assertRegex(rep, f'{self.tas.int_seq2.uuid}')

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 2
        )

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
        self.assertRegex(
            rep, 'ProvDAG representing the provenance.*Artifacts.*'
        )
        self.assertRegex(rep, f'{self.tas.single_int.uuid}')
        self.assertRegex(rep, f'{self.tas.int_seq1.uuid}')
        self.assertRegex(rep, f'{self.tas.mapping1.uuid}')

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 3
        )

    def test_union_self_missing_checksums_md5(self):
        unioned_dag = ProvDAG.union(
            [self.tas.dag_missing_md5, self.tas.single_int.dag]
        )

        self.assertRegex(
            repr(unioned_dag),
            'ProvDAG representing the provenance.*Artifacts.*'
            f'{self.tas.single_int.uuid}'
        )

        # The ChecksumDiff==None from the tinkered dag gets ignored...
        self.assertEqual(unioned_dag.checksum_diff, ChecksumDiff({}, {}, {}))

        # ...but this should make clear that the provenance is bad
        # (or that the user opted out of validation)
        self.assertEqual(
            unioned_dag.provenance_is_valid, ValidationCode.INVALID)

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1
        )

    def test_union_other_missing_checksums_md5(self):
        '''
        Tests unions of v5 dags where the other ProvDAG is missing its
        checksums.md5 but the calling ProvDAG is not
        '''
        unioned_dag = ProvDAG.union([self.tas.single_int.dag,
                                     self.tas.dag_missing_md5])

        self.assertRegex(repr(unioned_dag),
                         'ProvDAG representing the provenance.*Artifacts.*'
                         f'{self.tas.single_int.uuid}')

        self.assertEqual(unioned_dag.checksum_diff, ChecksumDiff({}, {}, {}))

        self.assertEqual(
            unioned_dag.provenance_is_valid, ValidationCode.INVALID
        )

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1
        )

    def test_union_both_missing_checksums_md5(self):
        '''
        Tests unions of v5 dags where both artifacts are missing their
        checksums.md5 files.
        '''
        unioned_dag = ProvDAG.union(
            [self.tas.dag_missing_md5, self.tas.dag_missing_md5])

        self.assertRegex(
            repr(unioned_dag),
            'ProvDAG representing the provenance.*Artifacts.*'
            f'{self.tas.single_int.uuid}'
        )

        # Both DAGs have NoneType checksum_diffs, so the ChecksumDiff==None
        self.assertEqual(unioned_dag.checksum_diff, None)

        self.assertEqual(
            unioned_dag.provenance_is_valid, ValidationCode.INVALID
        )

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1
        )

    def test_union_v0_v1_archives(self):
        unioned_dag = ProvDAG.union(
            [self.tas.table_v0.dag, self.tas.concated_ints_v2.dag]
        )

        self.assertIn(f'{self.tas.table_v0.uuid}', repr(unioned_dag))
        self.assertIn(f'{self.tas.concated_ints_v2.uuid}', repr(unioned_dag))

        self.assertEqual(
            unioned_dag.provenance_is_valid, ValidationCode.PREDATES_CHECKSUMS
        )

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 2
        )

        self.assertFalse(
            unioned_dag.node_has_provenance(self.tas.table_v0.uuid)
        )
        self.assertTrue(
            unioned_dag.node_has_provenance(self.tas.concated_ints_v2.uuid)
        )

    def test_union_v3_v5_archives(self):
        unioned_dag = ProvDAG.union(
            [self.tas.concated_ints_v3.dag, self.tas.concated_ints_v5.dag]
        )

        self.assertIn(f'{self.tas.concated_ints_v3.uuid}', repr(unioned_dag))
        self.assertIn(f'{self.tas.concated_ints_v5.uuid}', repr(unioned_dag))

        self.assertEqual(
            unioned_dag.provenance_is_valid,
            ValidationCode.PREDATES_CHECKSUMS
        )

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 2
        )

        self.assertTrue(
            unioned_dag.node_has_provenance(self.tas.concated_ints_v3.uuid)
        )
        self.assertTrue(
            unioned_dag.node_has_provenance(self.tas.concated_ints_v5.uuid)
        )

    def test_union_v5_v6_archives(self):
        unioned_dag = ProvDAG.union(
            [self.tas.concated_ints_v5.dag, self.tas.concated_ints_v6.dag]
        )

        self.assertIn(f'{self.tas.concated_ints_v5.uuid}', repr(unioned_dag))
        self.assertIn(f'{self.tas.concated_ints_v6.uuid}', repr(unioned_dag))

        self.assertEqual(
            unioned_dag.provenance_is_valid, ValidationCode.VALID
        )

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 2
        )

        self.assertTrue(
            unioned_dag.node_has_provenance(self.tas.concated_ints_v5.uuid)
        )
        self.assertTrue(
            unioned_dag.node_has_provenance(self.tas.concated_ints_v6.uuid)
        )

    def test_dag_is_superset(self):
        '''
        Tests union of three dags, where one dag is a proper superset of the
        others. We expect three _parsed_artifact_uuids, one terminal uuid,
        and one weakly_connected_component.
        '''
        unioned_dag = ProvDAG.union([
            self.tas.concated_ints.dag,
            self.tas.int_seq1.dag,
            self.tas.int_seq2.dag
        ])

        self.assertIn(
            self.tas.int_seq1.uuid, unioned_dag._parsed_artifact_uuids
        )
        self.assertIn(
            self.tas.int_seq2.uuid, unioned_dag._parsed_artifact_uuids
        )
        self.assertIn(
            self.tas.concated_ints.uuid, unioned_dag._parsed_artifact_uuids
        )
        self.assertEqual(len(unioned_dag._parsed_artifact_uuids), 3)

        self.assertEqual(len(unioned_dag.terminal_uuids), 1)
        self.assertEqual(
            unioned_dag.terminal_uuids, {self.tas.concated_ints.uuid}
        )

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1
        )

        # == tests identity of objects in memory, so we need is_isomorphic
        self.assertTrue(
            nx.is_isomorphic(self.tas.concated_ints.dag.dag, unioned_dag.dag)
        )
        self.assertEqual(self.tas.concated_ints.dag, unioned_dag)

    def test_three_artifacts_two_terminal_uuids(self):
        '''
        Tests union of four dags, where two artifacts have shared parents but
        no direct relationship with eachother. We expect four
        _parsed_artifact_uuids, two terminal uuids, and one
        weakly_connected_component.
        '''
        unioned_dag = ProvDAG.union([
            self.tas.int_seq1.dag,
            self.tas.int_seq2.dag,
            self.tas.concated_ints.dag,
            self.tas.other_concated_ints.dag
        ])

        self.assertIn(
            self.tas.int_seq1.uuid, unioned_dag._parsed_artifact_uuids
        )
        self.assertIn(
            self.tas.int_seq2.uuid, unioned_dag._parsed_artifact_uuids
        )
        self.assertIn(
            self.tas.concated_ints.uuid, unioned_dag._parsed_artifact_uuids
        )
        self.assertIn(
            self.tas.other_concated_ints.uuid,
            unioned_dag._parsed_artifact_uuids
        )
        self.assertEqual(len(unioned_dag._parsed_artifact_uuids), 4)

        self.assertEqual(len(unioned_dag.terminal_uuids), 2)
        self.assertEqual(
            unioned_dag.terminal_uuids,
            {self.tas.concated_ints.uuid, self.tas.other_concated_ints.uuid}
        )

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1
        )

    def test_one_analysis_two_artifacts(self):
        '''
        In this set of test archives, both artifacts are derived from the same
        parents, so should produce one connected DAG even though we are
        missing the parent artifacts used to create them
        '''
        unioned_dag = ProvDAG.union([
            self.tas.concated_ints.dag,
            self.tas.other_concated_ints.dag
        ])

        self.assertIn(
            self.tas.concated_ints.uuid, unioned_dag._parsed_artifact_uuids
        )
        self.assertIn(
            self.tas.other_concated_ints.uuid,
            unioned_dag._parsed_artifact_uuids
        )
        self.assertEqual(len(unioned_dag._parsed_artifact_uuids), 2)

        self.assertEqual(len(unioned_dag.terminal_uuids), 2)
        self.assertEqual(
            unioned_dag.terminal_uuids,
            {self.tas.concated_ints.uuid, self.tas.other_concated_ints.uuid}
        )

        self.assertEqual(
            nx.number_weakly_connected_components(unioned_dag.dag), 1
        )

    def test_no_checksum_validation_on_intact_artifact(self):
        no_validation_dag = ProvDAG(
            self.tas.concated_ints.filepath, validate_checksums=False
        )

        self.assertEqual(len(no_validation_dag.terminal_uuids), 1)
        terminal_uuid, *_ = no_validation_dag.terminal_uuids
        self.assertEqual(terminal_uuid, self.tas.concated_ints.uuid)

        self.assertEqual(len(no_validation_dag), 3)
        self.assertEqual(
            no_validation_dag.provenance_is_valid,
            ValidationCode.VALIDATION_OPTOUT
        )
        self.assertEqual(no_validation_dag.checksum_diff, None)

    def test_no_checksum_validation_missing_checksums_md5(self):
        with generate_archive_with_file_removed(
            self.tas.concated_ints.filepath,
            self.tas.concated_ints.uuid,
            'checksums.md5'
        ) as altered_archive:
            dag = ProvDAG(altered_archive, validate_checksums=False)

            self.assertEqual(
                dag.provenance_is_valid, ValidationCode.VALIDATION_OPTOUT
            )

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

                    expected = (f'(?s)Malformed.*{file}.*{uuid}.*corrupt.*')
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
            TypeError, f"EmptyParser.*{fn} is not None"
        ):
            EmptyParser.get_parser(fp)

    def test_parse_a_nonetype(self):
        '''
        tests that we can actually create empty ProvDAGs
        '''
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
        cls.tas = TestArtifacts()

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_get_parser(self):
        for archive in self.tas.all_artifact_versions:
            parser = ProvDAGParser.get_parser(archive.dag)
            self.assertIsInstance(parser, ProvDAGParser)

    def test_get_parser_input_data_not_a_provdag(self):
        fn = 'not_a_zip.txt'
        fp = os.path.join(self.tas.tempdir, fn)
        with self.assertRaisesRegex(
                TypeError, f"ProvDAGParser.*{fn} is not a ProvDAG"):
            ProvDAGParser.get_parser(fp)

    def test_parse_a_provdag(self):
        parser = ProvDAGParser()
        for archive in self.tas.all_artifact_versions:
            dag = archive.dag
            parsed = parser.parse_prov(Config(), dag)
            self.assertIsInstance(parsed, ParserResults)
            self.assertEqual(
                parsed.parsed_artifact_uuids, dag._parsed_artifact_uuids
            )
            self.assertTrue(nx.is_isomorphic(parsed.prov_digraph, dag.dag))
            self.assertEqual(
                parsed.provenance_is_valid, dag.provenance_is_valid
            )
            self.assertEqual(parsed.checksum_diff, dag.checksum_diff)


class SelectParserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_correct_parser_type(self):
        empty = select_parser(None)
        self.assertIsInstance(empty, EmptyParser)

        archive = select_parser(self.tas.concated_ints.filepath)
        self.assertIsInstance(archive, ArchiveParser)

        dag = ProvDAG()
        pdag = select_parser(dag)
        self.assertIsInstance(pdag, ProvDAGParser)

        # check dir_fp as fp
        test_dir = os.path.join(self.tempdir, 'parse_dir_test')
        os.mkdir(test_dir)
        dir_fp = pathlib.Path(test_dir)
        dir_p = select_parser(dir_fp)
        self.assertIsInstance(dir_p, DirectoryParser)

        # check dir_fp as str
        dir_fp_str = str(dir_fp)
        dir_p = select_parser(dir_fp_str)
        self.assertIsInstance(dir_p, DirectoryParser)

    def test_correct_archive_parser_version(self):
        parsers = [
            ParserV0, ParserV1, ParserV2, ParserV3, ParserV4, ParserV5,
            ParserV6
        ]
        for archive, parser in zip(self.tas.all_artifact_versions, parsers):
            handler = select_parser(archive.filepath)
            self.assertEqual(type(handler), parser)


class ParseProvenanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir
        cls.cfg = Config()

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_parse_with_artifact_parser(self):
        uuid = self.tas.concated_ints.uuid
        fp = self.tas.concated_ints.filepath
        parser_results = parse_provenance(self.cfg, fp)

        self.assertIsInstance(parser_results, ParserResults)
        p_a_uuids = parser_results.parsed_artifact_uuids
        self.assertIsInstance(p_a_uuids, set)
        self.assertIsInstance(next(iter(p_a_uuids)), str)
        self.assertEqual(len(parser_results.prov_digraph), 3)
        self.assertIn(uuid, parser_results.prov_digraph)
        self.assertIsInstance(
            parser_results.prov_digraph.nodes[uuid]['node_data'],
            ProvNode
        )
        self.assertEqual(
            parser_results.provenance_is_valid, ValidationCode.VALID
        )
        self.assertEqual(
            parser_results.checksum_diff, ChecksumDiff({}, {}, {})
        )

    def test_parse_with_provdag_parser(self):
        uuid = self.tas.concated_ints.uuid
        dag = self.tas.concated_ints.dag
        parser_results = parse_provenance(self.cfg, dag)

        self.assertIsInstance(parser_results, ParserResults)
        p_a_uuids = parser_results.parsed_artifact_uuids
        self.assertIsInstance(p_a_uuids, set)
        self.assertIsInstance(next(iter(p_a_uuids)), str)
        self.assertEqual(len(parser_results.prov_digraph), 3)
        self.assertIn(uuid, parser_results.prov_digraph)
        self.assertIsInstance(
            parser_results.prov_digraph.nodes[uuid]['node_data'],
            ProvNode
        )
        self.assertEqual(
            parser_results.provenance_is_valid, ValidationCode.VALID
        )
        self.assertEqual(
            parser_results.checksum_diff, ChecksumDiff({}, {}, {})
        )

    def test_parse_with_empty_parser(self):
        res = parse_provenance(self.cfg, None)
        self.assertIsInstance(res, ParserResults)
        self.assertEqual(res.parsed_artifact_uuids, set())
        self.assertTrue(nx.is_isomorphic(res.prov_digraph, nx.DiGraph()))
        self.assertEqual(res.provenance_is_valid, ValidationCode.VALID)
        self.assertEqual(res.checksum_diff, None)

    def test_parse_with_directory_parser(self):
        # Non-recursive
        parse_dir_fp = os.path.join(self.tempdir, 'parse_dir')
        os.mkdir(parse_dir_fp)
        concated_ints_path = os.path.join(parse_dir_fp, 'concated-ints.qza')
        shutil.copy(self.tas.concated_ints.filepath, concated_ints_path)

        res = parse_provenance(self.cfg, parse_dir_fp)
        self.assertEqual(self.cfg.recurse, False)
        self.assertIsInstance(res, ParserResults)
        concated_ints_uuid = self.tas.concated_ints.uuid
        self.assertEqual(res.parsed_artifact_uuids, {concated_ints_uuid})
        self.assertEqual(len(res.prov_digraph), 3)
        self.assertEqual(res.provenance_is_valid, ValidationCode.VALID)
        self.assertEqual(res.checksum_diff, ChecksumDiff({}, {}, {}))

        # Recursive
        inner_dir_path = os.path.join(parse_dir_fp, 'inner-dir')
        os.mkdir(inner_dir_path)
        mapping_path = os.path.join(inner_dir_path, 'mapping1.qza')
        shutil.copy(self.tas.mapping1.filepath, mapping_path)

        self.cfg.recurse = True
        self.assertEqual(self.cfg.recurse, True)
        res = parse_provenance(self.cfg, parse_dir_fp)
        self.assertIsInstance(res, ParserResults)
        mapping_uuid = self.tas.mapping1.uuid
        self.assertEqual(
            res.parsed_artifact_uuids, {concated_ints_uuid, mapping_uuid}
        )
        self.assertEqual(len(res.prov_digraph), 4)
        self.assertEqual(res.provenance_is_valid, ValidationCode.VALID)
        self.assertEqual(res.checksum_diff, ChecksumDiff({}, {}, {}))

    def test_parse_with_directory_parser_bad_dir_path(self):
        dir_fp = os.path.join(self.tempdir, 'fake_dir')
        with self.assertRaisesRegex(Exception, 'not a valid dir'):
            parse_provenance(self.cfg, dir_fp)

    def test_no_correct_parser_found_error(self):
        input_data = {'this': 'is not parseable'}
        with self.assertRaisesRegex(
            UnparseableDataError,
            f'(?s)Input data {input_data}.*not supported.*'
            'ArchiveParser expects a string or pathlib.PosixPath.*'
            'DirectoryParser.*expects a directory.*'
            'ProvDAGParser.*is not a ProvDAG.*'
            'EmptyParser.*is not None'
        ):
            select_parser(input_data)


class DirectoryParserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir
        cls.cfg = Config()

        parse_dir_fp = os.path.join(cls.tempdir, 'parse-dir')
        os.mkdir(parse_dir_fp)
        concated_ints_path = os.path.join(parse_dir_fp, 'concated-ints.qza')
        shutil.copy(cls.tas.concated_ints.filepath, concated_ints_path)
        int_seq_path = os.path.join(parse_dir_fp, 'int-seq1.qza')
        shutil.copy(cls.tas.int_seq1.filepath, int_seq_path)

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_parse_empty_dir(self):
        empty_dir_path = os.path.join(self.tempdir, 'empty-dir')
        os.mkdir(empty_dir_path)

        with self.assertRaisesRegex(
            ValueError, f'No .qza or .qzv files.*{empty_dir_path}'
        ):
            ProvDAG(empty_dir_path)

    def test_directory_parser_works_regardless_trailing_slash(self):
        dag = ProvDAG(self.tempdir + '/parse-dir')
        dag2 = ProvDAG(self.tempdir + '/parse-dir/')
        self.assertEqual(dag, dag2)

    def test_directory_parser_captures_all_parsed_artifact_uuids(self):
        '''
        The test dir contains a concatenated ints artifact and an int sequence.
        The concatenated ints artifact is a true subset of the int sequence,
        and is not re-parsed as implemented.

        The resulting dag should be aware of both user-passed artifacts,
        regardless of whether we actually parse the int sequence.
        '''
        dag = ProvDAG(self.tempdir + '/parse-dir/', recurse=True)
        self.assertEqual(
            dag._parsed_artifact_uuids,
            {self.tas.concated_ints.uuid, self.tas.int_seq1.uuid}
        )

    def test_directory_parser_idempotent_with_parse_and_union(self):
        # Non-recursive
        concated_ints_path = os.path.join(
            self.tempdir, 'parse-dir', 'concated-ints.qza'
        )
        int_seq_path = os.path.join(self.tempdir, 'parse-dir', 'int-seq1.qza')
        concated_ints_dag = ProvDAG(concated_ints_path)
        int_seq_dag = ProvDAG(int_seq_path)
        union_dag = ProvDAG.union([concated_ints_dag, int_seq_dag])
        parse_dir_path = os.path.join(self.tempdir, 'parse-dir')
        dir_dag = ProvDAG(parse_dir_path)
        self.assertEqual(union_dag, dir_dag)

        # Recursive
        inner_dir = os.path.join(self.tempdir, 'parse-dir', 'inner-dir')
        os.mkdir(inner_dir)
        mapping_path = os.path.join(inner_dir, 'mapping1.qza')
        shutil.copy(self.tas.mapping1.filepath, mapping_path)

        mapping_dag = ProvDAG(mapping_path)
        inner_dir_union_dag = ProvDAG.union([
            concated_ints_dag, int_seq_dag, mapping_dag
        ])
        recursive_dir_dag = ProvDAG(parse_dir_path, recurse=True)
        self.assertEqual(inner_dir_union_dag, recursive_dir_dag)

    def test_directory_parser_multiple_imports(self):
        outer_path = os.path.join(self.tempdir, 'mutliple-import-outer')
        os.mkdir(outer_path)
        inner_path = os.path.join(outer_path, 'inner')
        os.mkdir(inner_path)

        shutil.copy(
            self.tas.mapping1.filepath,
            os.path.join(outer_path, 'mapping1.qza')
        )
        shutil.copy(
            self.tas.mapping2.filepath,
            os.path.join(outer_path, 'mapping2.qza')
        )
        shutil.copy(
            self.tas.mapping1.filepath,
            os.path.join(inner_path, 'duplicate-mapping1.qza')
        )
        shutil.copy(
            self.tas.mapping2.filepath,
            os.path.join(inner_path, 'duplicate-mapping2.qza')
        )

        inner_dag = ProvDAG(inner_path)
        self.assertEqual(len(inner_dag), 2)
        self.assertIn(self.tas.mapping1.uuid, inner_dag.dag)
        self.assertIn(self.tas.mapping2.uuid, inner_dag.dag)

        outer_dag = ProvDAG(outer_path)
        self.assertEqual(len(inner_dag), 2)
        self.assertIn(self.tas.mapping1.uuid, outer_dag.dag)
        self.assertIn(self.tas.mapping2.uuid, outer_dag.dag)
        self.assertEqual(inner_dag, outer_dag)

        recursive_outer_dag = ProvDAG(outer_path, recurse=True)
        self.assertEqual(len(inner_dag), 2)
        self.assertIn(self.tas.mapping1.uuid, recursive_outer_dag.dag)
        self.assertIn(self.tas.mapping2.uuid, recursive_outer_dag.dag)
        self.assertEqual(inner_dag, outer_dag, recursive_outer_dag)

    def test_verbose(self):
        buffer = io.StringIO()
        concated_ints = 'parse-dir/concated-ints.qza'
        int_seq = 'parse-dir/int-seq1.qza'
        with redirect_stdout(buffer):
            dag = ProvDAG(
                os.path.join(self.tempdir, 'parse-dir'),
                verbose=True,
                recurse=True
            )

        self.assertEqual(dag.cfg.verbose, True)

        stdout_log = buffer.getvalue()
        self.assertRegex(stdout_log, f'parsing.*{concated_ints}')
        self.assertRegex(stdout_log, f'parsing.*{int_seq}')
