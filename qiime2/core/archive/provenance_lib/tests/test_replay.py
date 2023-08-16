import bibtexparser as bp
import networkx as nx
import os
import pandas as pd
import pathlib
import re
import shutil
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from qiime2 import Artifact
from qiime2.sdk import PluginManager
from qiime2.sdk.usage import Usage, UsageVariable
from q2cli.core.usage import CLIUsageVariable
from qiime2.plugins import ArtifactAPIUsageVariable

from ..parse import ProvDAG
from ..replay import (
    ActionCollections, BibContent, NamespaceCollections, ReplayConfig,
    UsageVarsDict,
    build_no_provenance_node_usage, build_import_usage, build_action_usage,
    build_usage_examples, collect_citations, dedupe_citations,
    dump_recorded_md_file, group_by_action, init_md_from_artifacts,
    init_md_from_md_file, init_md_from_recorded_md, param_is_metadata_column,
    replay_provenance, uniquify_action_name, replay_citations,
    replay_supplement,
    SUPPORTED_USAGE_DRIVERS,
    )
from .testing_utilities import (
    CustomAssertions, DATA_DIR, TEST_DATA, TestArtifacts
)
from ..util import camel_to_snake
from ...provenance import MetadataInfo


class UsageVarsDictTests(unittest.TestCase):
    def test_uniquify(self):
        collision_val = 'emp_single_end_sequences'
        unique_val = 'some_prime'
        ns = UsageVarsDict({'123': collision_val})
        self.assertEqual(ns.data, {'123': 'emp_single_end_sequences_0'})
        ns.update({'456': collision_val, 'unique': unique_val})
        self.assertEqual(ns['456'], 'emp_single_end_sequences_1')
        self.assertEqual(ns['unique'], 'some_prime_0')
        ns['789'] = collision_val
        self.assertEqual(ns.pop('789'), 'emp_single_end_sequences_2')

    def test_get_key(self):
        ns = UsageVarsDict({'123': 'some_name'})
        self.assertEqual('123', ns.get_key('some_name_0'))
        with self.assertRaisesRegex(
            KeyError, "passed value 'fake_key' does not exist"
        ):
            ns.get_key('fake_key')


class NamespaceCollectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_add_usage_var_workflow(self):
        """
        Smoke tests a common workflow with this data structure
        - Create a unique variable name by adding to .usg_var_namespace
        - Create a UsageVariable with that name
        - use the name to get the UUID (when we have Results, we have no UUIDs)
        - add the correctly-named UsageVariable to .usg_vars
        """
        use = Usage()
        uuid = self.tas.concated_ints.uuid
        base_name = 'concated_ints'
        exp_name = base_name + '_0'
        ns = NamespaceCollections()
        ns.usg_var_namespace.update({uuid: base_name})
        self.assertEqual(ns.usg_var_namespace[uuid], exp_name)

        def factory():  # pragma: no cover
            return Artifact.load(self.tas.concated_ints.filepath)
        u_var = use.init_artifact(ns.usg_var_namespace[uuid], factory)
        self.assertEqual(u_var.name, exp_name)

        actual_uuid = ns.usg_var_namespace.get_key(u_var.name)
        self.assertEqual(actual_uuid, uuid)

        ns.usg_vars[uuid] = u_var
        self.assertIsInstance(ns.usg_vars[uuid], UsageVariable)
        self.assertEqual(ns.usg_vars[uuid].name, exp_name)


class ReplayProvenanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_replay_from_fp(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'rendered.txt'
            out_fn = str(out_fp)
            in_fn = self.tas.concated_ints_with_md.filepath
            replay_provenance(in_fn, out_fn, 'python3')

            self.assertTrue(out_fp.is_file())

            with open(out_fn, 'r') as fp:
                rendered = fp.read()

            self.assertIn('from qiime2 import Artifact', rendered)
            self.assertIn('from qiime2 import Metadata', rendered)
            self.assertIn(
                'import qiime2.plugins.dummy_plugin.actions as '
                'dummy_plugin_actions',
                rendered
            )
            self.assertIn('mapping_0 = Artifact.import_data(', rendered)
            self.assertRegex(rendered,
                             'The following command.*additional metadata')
            self.assertIn('mapping_0.view(Metadata)', rendered)
            self.assertIn('dummy_plugin_actions.identity_with_metadata',
                          rendered)
            self.assertIn('dummy_plugin_actions.concatenate_ints', rendered)

    def test_replay_from_fp_use_md_without_parse(self):
        in_fp = self.tas.concated_ints.filepath
        with self.assertRaisesRegex(
                ValueError, "Metadata not parsed for replay. Re-run"
        ):
            replay_provenance(
                in_fp, 'unused_fp', 'python3',
                parse_metadata=False, use_recorded_metadata=True
            )

    def test_replay_dump_md_without_parse(self):
        in_fp = self.tas.concated_ints.filepath
        with self.assertRaisesRegex(
                ValueError, "(?s)Metadata not parsed,.*dump_recorded_meta"
        ):
            replay_provenance(
                in_fp, 'unused_fp', 'python3',
                parse_metadata=False, dump_recorded_metadata=True
            )

    def test_replay_md_out_fp_without_parse(self):
        in_fp = self.tas.concated_ints.filepath
        with self.assertRaisesRegex(
                ValueError, "(?s)Metadata not parsed,.*not.*metadata output"
        ):
            replay_provenance(
                in_fp, 'unused_fp', 'python3', parse_metadata=False,
                dump_recorded_metadata=False,
                md_out_fp='/user/dumb/some_filepath'
            )

    def test_replay_use_md_without_dump_md(self):
        in_fp = self.tas.concated_ints.filepath
        with self.assertRaisesRegex(
                NotImplementedError,
                "(?s)uses.*metadata.*must.*written to disk"
        ):
            replay_provenance(
                in_fp, 'unused_fp', 'python3', use_recorded_metadata=True,
                dump_recorded_metadata=False
            )

    def test_replay_from_provdag(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'rendered.txt'
            out_fn = str(out_fp)
            dag = self.tas.concated_ints_with_md.dag
            replay_provenance(dag, out_fn, 'python3')

            self.assertTrue(out_fp.is_file())

            with open(out_fn, 'r') as fp:
                rendered = fp.read()

            self.assertIn('from qiime2 import Artifact', rendered)
            self.assertIn('from qiime2 import Metadata', rendered)
            self.assertIn(
                'import qiime2.plugins.dummy_plugin.actions as '
                'dummy_plugin_actions',
                rendered
            )
            self.assertIn('mapping_0 = Artifact.import_data(', rendered)
            self.assertRegex(rendered,
                             'The following command.*additional metadata')
            self.assertIn('mapping_0.view(Metadata)', rendered)
            self.assertIn('dummy_plugin_actions.identity_with_metadata',
                          rendered)
            self.assertIn('dummy_plugin_actions.concatenate_ints', rendered)

    def test_replay_from_provdag_use_md_without_parse(self):
        dag = ProvDAG(self.tas.concated_ints_with_md.filepath,
                      validate_checksums=False,
                      parse_metadata=False)
        with self.assertRaisesRegex(
                ValueError, "Metadata not parsed for replay"
        ):
            replay_provenance(dag, 'unused', 'python3',
                              use_recorded_metadata=True)

    def test_replay_from_provdag_ns_collision(self):
        """
        This artifact's dag contains a few results with the output-name
        filtered-table, so is a good check for namespace collisions if
        we're not uniquifying variable names properly.
        """
        with tempfile.TemporaryDirectory() as tempdir:
            self.tas.concated_ints.artifact.save(
                os.path.join(tempdir, 'c1.qza'))
            self.tas.other_concated_ints.artifact.save(
                os.path.join(tempdir, 'c2.qza'))
            dag = ProvDAG(tempdir)

        drivers = ['python3', 'cli']
        exp = ['concatenated_ints_0', 'concatenated_ints_1']
        for driver in drivers:
            with tempfile.TemporaryDirectory() as tempdir:
                out_path = pathlib.Path(tempdir) / 'ns_coll.txt'
                replay_provenance(dag, out_path, driver)

                with open(out_path, 'r') as fp:
                    rendered = fp.read()
                    for name in exp:
                        if driver == 'cli':
                            name = name.replace('_', '-')
                        self.assertIn(name, rendered)

    def test_replay_optional_param_is_none(self):
        dag = self.tas.int_seq_optional_input.dag
        with tempfile.TemporaryDirectory() as tempdir:
            out_path = pathlib.Path(tempdir) / 'ns_coll.txt'

            replay_provenance(dag, out_path, 'python3')
            with open(out_path, 'r') as fp:
                rendered = fp.read()
            self.assertIn('ints=int_sequence1_0', rendered)
            self.assertIn('num1=', rendered)
            self.assertNotIn('optional1=', rendered)
            self.assertNotIn('num2=', rendered)

            replay_provenance(dag, out_path, 'cli')
            with open(out_path, 'r') as fp:
                rendered = fp.read()
            self.assertIn('--i-ints int-sequence1-0.qza', rendered)
            self.assertIn('--p-num1', rendered)
            self.assertNotIn('--i-optional1', rendered)
            self.assertNotIn('--p-num2', rendered)

    # TODO: not working after removing `action_patch`, need to remember why
    # we thought that wasn't an issue
    def test_replay_untracked_output_names(self):
        """
        In this artifact, some nodes don't track output names. As
        a result, replay could fail when Usage.action tries to look up the
        output-name for those results in the plugin manager's record of the
        actual action's signature. We've monkeypatched Usage.action here to
        allow replay to proceed.
        """
        dag = self.tas.concated_ints_no_output_names.dag
        # If we rendered the final action correctly, then nothing blew up.
        with tempfile.TemporaryDirectory() as tempdir:
            out_path = pathlib.Path(tempdir) / 'ns_coll.txt'

            replay_provenance(dag, out_path, 'python3')
            with open(out_path, 'r') as fp:
                rendered = fp.read()
            self.assertIn('dummy_plugin_actions.concatenate_ints', rendered)

            replay_provenance(dag, out_path, 'cli')
            with open(out_path, 'r') as fp:
                rendered = fp.read()
            self.assertIn('concatenate_ints', rendered)
            self.assertIn('--o-concatenated-ints', rendered)

        # dag = ProvDAG(os.path.join(DATA_DIR, 'heatmap.qzv'))
        # drivers = ['python3', 'cli']
        # # If we rendered the final action correctly, then nothing blew up.
        # exp = {
        #     'python3':
        #         '(?s)action_results.*classifier_actions.classify_samples.*'
        #         'heatmap_0_viz = action_results.heatmap',
        #     'cli': '(?s)qiime sample-classifier classify-samples.*'
        #            '--o-heatmap heatmap-0.qzv'}
        # for driver in drivers:
        #     with tempfile.TemporaryDirectory() as tmpdir:
        #         out_path = pathlib.Path(tmpdir) / 'ns_coll.txt'
        #         replay_provenance(dag, out_path, driver)

        #         with open(out_path, 'r') as fp:
        #             rendered = fp.read()
        #     self.assertRegex(rendered, exp[driver])


class ReplayProvDAGDirectoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_directory_replay_multiple_imports(self):
        """
        The directory being parsed here contains two pairs of duplicates,
        and should replay as only two import statements.
        """
        outer_dir = os.path.join(self.tempdir, 'outer')
        inner_dir = os.path.join(outer_dir, 'inner')
        os.makedirs(inner_dir)

        for artifact in self.tas.single_int, self.tas.single_int2:
            for dir_ in inner_dir, outer_dir:
                artifact.artifact.save(
                    os.path.join(dir_, f'{artifact.name}.qza')
                )

        dir_dag = ProvDAG(outer_dir)
        self.assertEqual(len(dir_dag._parsed_artifact_uuids), 2)
        self.assertIn(self.tas.single_int.uuid, dir_dag.dag)
        self.assertIn(self.tas.single_int2.uuid, dir_dag.dag)

        exp_1 = (
            '(?s)from qiime2 import Artifact.*'
            'single_int_0 = Artifact.import_data.*'
            '<your data here>.*'
        )
        exp_2 = (
            '(?s)from qiime2 import Artifact.*'
            'single_int_1 = Artifact.import_data.*'
            '<your data here>.*'
        )

        with tempfile.TemporaryDirectory() as tempdir:
            out_path = pathlib.Path(tempdir) / 'rendered.txt'
            replay_provenance(dir_dag, out_path, 'python3')
            self.assertTrue(out_path.is_file())

            with open(out_path, 'r') as fp:
                rendered = fp.read()
                self.assertRegex(rendered, exp_1)
                self.assertRegex(rendered, exp_2)


class BuildUsageExamplesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir
        cls.pm = PluginManager()

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    @patch('qiime2.core.archive.provenance_lib.replay.build_action_usage')
    @patch('qiime2.core.archive.provenance_lib.replay.build_import_usage')
    @patch('qiime2.core.archive.provenance_lib.replay.'
           'build_no_provenance_node_usage')
    def test_build_usage_examples(self, n_p_builder, imp_builder, act_builder):
        dag = self.tas.concated_ints_with_md.dag
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=self.pm)
        build_usage_examples(dag, cfg)

        n_p_builder.assert_not_called()
        self.assertEqual(imp_builder.call_count, 3)
        self.assertEqual(act_builder.call_count, 2)

    @patch('qiime2.core.archive.provenance_lib.replay.build_action_usage')
    @patch('qiime2.core.archive.provenance_lib.replay.build_import_usage')
    @patch('qiime2.core.archive.provenance_lib.replay.'
           'build_no_provenance_node_usage')
    def test_build_usage_examples_lone_v0(
            self, n_p_builder, imp_builder, act_builder
    ):
        uuid = self.tas.table_v0.uuid
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{uuid}.*prior.*incomplete'
        ):
            dag = ProvDAG(self.tas.table_v0.filepath)

        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=self.pm)
        build_usage_examples(dag, cfg)

        # This is a single v0 archive, so should have only one np node
        n_p_builder.assert_called_once()
        imp_builder.assert_not_called()
        act_builder.assert_not_called()

    @patch('qiime2.core.archive.provenance_lib.replay.build_action_usage')
    @patch('qiime2.core.archive.provenance_lib.replay.build_import_usage')
    @patch('qiime2.core.archive.provenance_lib.replay.'
           'build_no_provenance_node_usage')
    def test_build_usage_examples_mixed(
            self, n_p_builder, imp_builder, act_builder
    ):
        mixed_dir = os.path.join(self.tempdir, 'mixed-dir')
        os.mkdir(mixed_dir)
        shutil.copy(self.tas.table_v0.filepath, mixed_dir)
        shutil.copy(self.tas.concated_ints_v6.filepath, mixed_dir)

        v0_uuid = self.tas.table_v0.uuid
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'
        ):
            dag = ProvDAG(mixed_dir)

        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=self.pm)
        build_usage_examples(dag, cfg)

        n_p_builder.assert_called_once()
        self.assertEqual(imp_builder.call_count, 2)
        act_builder.assert_called_once()

    @patch('qiime2.core.archive.provenance_lib.replay.build_action_usage')
    @patch('qiime2.core.archive.provenance_lib.replay.build_import_usage')
    @patch('qiime2.core.archive.provenance_lib.replay.'
           'build_no_provenance_node_usage')
    def test_build_usage_examples_big(
            self, n_p_builder, imp_builder, act_builder):

        many_dir = os.path.join(self.tempdir, 'many-dir')
        os.mkdir(many_dir)
        shutil.copy(self.tas.concated_ints_with_md.filepath, many_dir)
        shutil.copy(self.tas.splitted_ints.filepath, many_dir)
        shutil.copy(self.tas.pipeline_viz.filepath, many_dir)

        dag = ProvDAG(many_dir)
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=self.pm)
        build_usage_examples(dag, cfg)

        n_p_builder.assert_not_called()
        self.assertEqual(imp_builder.call_count, 3)
        self.assertEqual(act_builder.call_count, 4)

    # TODO: broken by Collections
    def test_build_action_usage_collection_of_inputs(self):
        """
        This artifact's root node is passed a collection of inputs. Collections
        need to be handled in a separate branch, so this exists for coverage.
        """
        dag = self.tas.int_from_collection.dag
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=self.pm)
        build_usage_examples(dag, cfg)
        list_line = 'tables=[table_0, table_1],'
        self.assertIn(list_line, cfg.use.render())

        # dag = ProvDAG(os.path.join(DATA_DIR, 'merged_tbls.qza'))
        # cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
        #                    use_recorded_metadata=False, pm=self.pm)
        # build_usage_examples(dag, cfg)
        # list_line = 'tables=[table_0, table_1],'
        # self.assertIn(list_line, cfg.use.render())


class MiscHelperFnTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir
        cls.pm = PluginManager()

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_uniquify_action_name(self):
        ns = set()
        p1 = 'dummy_plugin'
        a1 = 'action_jackson'
        p2 = 'dummy_plugin'
        a2 = 'missing_in_action'
        unique1 = uniquify_action_name(p1, a1, ns)
        self.assertEqual(unique1, 'dummy_plugin_action_jackson_0')
        unique2 = uniquify_action_name(p2, a2, ns)
        self.assertEqual(unique2, 'dummy_plugin_missing_in_action_0')
        duplicate = uniquify_action_name(p1, a1, ns)
        self.assertEqual(duplicate, 'dummy_plugin_action_jackson_1')

    def test_param_is_metadata_col(self):
        """
        Assumes q2-demux and q2-diversity are installed in the active env.
        """
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['cli'](),
                           use_recorded_metadata=False, pm=self.pm)

        actual = param_is_metadata_column(
            cfg, 'metadata', 'dummy_plugin', 'identity_with_metadata_column'
        )
        self.assertTrue(actual)

        actual = param_is_metadata_column(
            cfg, 'int1', 'dummy_plugin', 'concatenate_ints'
        )
        self.assertFalse(actual)

        with self.assertRaisesRegex(KeyError, "No action.*registered.*"):
            param_is_metadata_column(
                cfg, 'ints', 'dummy_plugin', 'young'
            )

        with self.assertRaisesRegex(KeyError, "No param.*registered.*"):
            param_is_metadata_column(
                cfg, 'thugger', 'dummy_plugin', 'split_ints'
            )

        with self.assertRaisesRegex(KeyError, "No plugin.*registered.*"):
            param_is_metadata_column(
                cfg, 'fake_param', 'dummy_hard', 'split_ints'
            )

    def test_dump_recorded_md_file_to_custom_dir(self):
        dag = self.tas.int_seq_with_md.dag
        uuid = self.tas.int_seq_with_md.uuid

        out_dir = 'custom_dir'
        provnode = dag.get_node_data(uuid)
        og_md = provnode.metadata['metadata']
        action_name = 'concatenate_ints_0'
        md_id = 'metadata'
        fn = 'metadata.tsv'

        with tempfile.TemporaryDirectory() as tempdir:
            cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['cli'](),
                               pm=self.pm, md_out_fp=(tempdir + '/' + out_dir))
            dump_recorded_md_file(cfg, provnode, action_name, md_id, fn)
            out_path = pathlib.Path(tempdir) / out_dir / action_name / fn

            self.assertTrue(out_path.is_file())

            dumped_df = pd.read_csv(out_path, sep='\t')
            pd.testing.assert_frame_equal(dumped_df, og_md)

            # If we run it again, it shouldn't overwrite 'recorded_metadata',
            # so we should have two files
            action_name_2 = 'concatenate_ints_1'
            md_id2 = 'metadata'
            fn2 = 'metadata_1.tsv'
            dump_recorded_md_file(cfg, provnode, action_name_2, md_id2, fn2)
            out_path2 = pathlib.Path(tempdir) / out_dir / action_name_2 / fn2

            # are both files where expected?
            self.assertTrue(out_path.is_file())
            self.assertTrue(out_path2.is_file())

    def test_dump_recorded_md_file_no_md(self):
        uuid = self.tas.table_v0.uuid
        dag = self.tas.table_v0.dag

        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           pm=self.pm)
        provnode = dag.get_node_data(uuid)
        action_name = 'old_action'
        md_id = 'metadata'
        fn = 'metadata.tsv'

        with self.assertRaisesRegex(
            ValueError, "should only be called.*if.*metadata"
        ):
            dump_recorded_md_file(cfg, provnode, action_name, md_id, fn)


class GroupByActionTests(unittest.TestCase):
    def test_g_b_a_w_provenance(self):
        self.maxDiff = None
        v5_dag = ProvDAG(TEST_DATA['5']['qzv_fp'])
        sorted_nodes = nx.topological_sort(v5_dag.collapsed_view)
        actual = group_by_action(v5_dag, sorted_nodes)
        exp = {
            '5cf3fd87-22ac-47ea-936d-576cc12f9110':
                {'a35830e1-4535-47c6-aa23-be295a57ee1c':
                    'emp_single_end_sequences', },
            '3d69c8d1-a1fa-4ab3-ac88-3a98da15b2d5':
                {'99fa3670-aa1a-45f6-ba8e-803c976a1163':
                    'per_sample_sequences', },
            'c0f74c1c-d596-4b1b-9aad-50101b4a9950':
                {'89af91c0-033d-4e30-8ac4-f29a3b407dc1': 'table',
                 '7ecf8954-e49a-4605-992e-99fcee397935':
                    'representative_sequences', },
            '4779bf7d-cae8-4ff2-a27d-4c97f581803f':
                {'bce3d09b-e296-4f2b-9af4-834db6412429': 'rooted_tree', },
            'aefba6e7-3dd1-45e5-8e6b-60e784062a5e':
                {'ffb7cee3-2f1f-4988-90cc-efd5184ef003':
                    'unweighted_unifrac_emperor', },
        }
        self.assertEqual(actual.std_actions, exp)
        self.assertEqual(actual.no_provenance_nodes, [])

    def test_g_b_a_no_provenance(self):
        # one v0 node
        v0_uuid = '0b8b47bd-f2f8-4029-923c-0e37a68340c3'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'):
            single_no_prov = ProvDAG(
                os.path.join(DATA_DIR, 'v0_uu_emperor.qzv'))
        sorted_nodes = nx.topological_sort(single_no_prov.collapsed_view)
        action_collections = group_by_action(single_no_prov, sorted_nodes)
        self.assertEqual(action_collections.std_actions, {})
        self.assertEqual(action_collections.no_provenance_nodes, [v0_uuid])

    def test_g_b_a_two_joined_no_prov_nodes(self):
        # Multiple no-prov nodes glued together into a single DAG
        v0_uuid = '0b8b47bd-f2f8-4029-923c-0e37a68340c3'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'):
            single_no_prov = ProvDAG(os.path.join(DATA_DIR,
                                     'v0_uu_emperor.qzv'))

        tbl_uuid = '89af91c0-033d-4e30-8ac4-f29a3b407dc1'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{tbl_uuid}.*prior.*incomplete'):
            v0_tbl = ProvDAG(os.path.join(DATA_DIR, 'v0_table.qza'))
        joined = ProvDAG.union([single_no_prov, v0_tbl])
        sorted_nodes = nx.topological_sort(joined.collapsed_view)
        action_collections = group_by_action(joined, sorted_nodes)
        self.assertEqual(action_collections.std_actions, {})
        self.assertEqual(action_collections.no_provenance_nodes,
                         [v0_uuid, tbl_uuid])

    def test_g_b_a_some_nodes_missing_provenance(self):
        # A dag parsed from a v1 archive with a v0 predecessor node
        act_id = 'c147dfbc-139a-4db0-ac17-b11948247f93'
        v1_uuid = '0b8b47bd-f2f8-4029-923c-0e37a68340c3'
        v0_uuid = '9f6a0f3e-22e6-4c39-8733-4e672919bbc7'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'):
            mixed = ProvDAG(os.path.join(DATA_DIR,
                            'mixed_v0_v1_uu_emperor.qzv'))
        sorted_nodes = nx.topological_sort(mixed.collapsed_view)
        action_collections = group_by_action(mixed, sorted_nodes)
        self.assertEqual(action_collections.std_actions,
                         {act_id: {v1_uuid: 'visualization'}})
        self.assertEqual(action_collections.no_provenance_nodes, [v0_uuid])


class InitializerTests(unittest.TestCase):
    def test_init_md_from_recorded_md(self):
        no_md_id = '89af91c0-033d-4e30-8ac4-f29a3b407dc1'
        has_md_id = '99fa3670-aa1a-45f6-ba8e-803c976a1163'
        dag = ProvDAG(os.path.join(DATA_DIR, 'v5_table.qza'))
        no_md_node = dag.get_node_data(no_md_id)
        md_node = dag.get_node_data(has_md_id)
        var_nm = 'per_sample_sequences_0_barcodes'
        param_nm = 'barcodes'
        # We expect the variable name has already been added to the namespace
        ns = UsageVarsDict({var_nm: param_nm})
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=pm)
        md_fn = 'demux_emp_single_0/barcodes_0'

        with self.assertRaisesRegex(ValueError, 'only.*call.*if.*metadata'):
            init_md_from_recorded_md(
                no_md_node, param_nm, var_nm, ns, cfg, md_fn)

        var = init_md_from_recorded_md(
            md_node, param_nm, var_nm, ns, cfg, md_fn)
        self.assertIsInstance(var, UsageVariable)
        self.assertEqual(var.var_type, 'column')

        rendered = cfg.use.render()
        self.assertRegex(rendered, 'from qiime2 import Metadata')
        exp = (r"barcodes_0_md = Metadata.load\('.*/"
               r"recorded_metadata/demux_emp_single_0/barcodes_0.tsv'")
        self.assertRegex(rendered, exp)

    def test_init_md_from_recorded_md_user_passed_fp(self):
        has_md_id = '99fa3670-aa1a-45f6-ba8e-803c976a1163'
        dag = ProvDAG(os.path.join(DATA_DIR, 'v5_table.qza'))
        md_node = dag.get_node_data(has_md_id)

        var_nm = 'per_sample_sequences_0_barcodes'
        param_nm = 'barcodes'
        # the variable name should already be added to the namespace
        ns = UsageVarsDict({var_nm: param_nm})
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=pm,
                           md_out_fp='./md_out')
        md_fn = 'test_a/test_o'

        var = init_md_from_recorded_md(
            md_node, param_nm, var_nm, ns, cfg, md_fn)
        self.assertIsInstance(var, UsageVariable)
        self.assertEqual(var.var_type, 'column')

        rendered = cfg.use.render()
        print(rendered)
        self.assertRegex(rendered, 'from qiime2 import Metadata')
        exp = (r"barcodes_0_md = Metadata.load\('.*md_out/test_a/test_o.tsv'")
        self.assertRegex(rendered, exp)

    def test_init_md_from_artifacts_no_artifacts(self):
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=pm)
        usg_vars = {}
        md_info = MetadataInfo([], relative_fp='hmm.tsv')
        with self.assertRaisesRegex(ValueError,
                                    "not.*used.*input_artifact_uuids.*empty"):
            init_md_from_artifacts(md_info, usg_vars, cfg)

    def test_init_md_from_artifacts_one_art(self):
        # This helper doesn't capture real data, so we're only smoke testing,
        # checking type, and confirming the repr looks reasonable.
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=pm)

        # We expect artifact vars have already been added to the namespace, so:
        a1 = cfg.use.init_artifact(name='thing1', factory=lambda: None)
        ns = NamespaceCollections(usg_vars={'uuid1': a1})

        md_info = MetadataInfo(['uuid1'], relative_fp='hmm.tsv')
        var = init_md_from_artifacts(md_info, ns, cfg)
        self.assertIsInstance(var, UsageVariable)
        self.assertEqual(var.var_type, 'metadata')
        rendered = var.use.render()
        exp = """from qiime2 import Metadata

thing1_a_0_md = thing1.view(Metadata)"""
        self.assertEqual(rendered, exp)

    def test_init_md_from_artifacts_many(self):
        # This helper doesn't capture real data, so we're only smoke testing,
        # checking type, and confirming the repr looks reasonable.
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=pm)

        # We expect artifact vars have already been added to the namespace, so:
        a1 = cfg.use.init_artifact(name='thing1', factory=lambda: None)
        a2 = cfg.use.init_artifact(name='thing2', factory=lambda: None)
        a3 = cfg.use.init_artifact(name='thing3', factory=lambda: None)
        ns = NamespaceCollections(
            usg_vars={'uuid1': a1, 'uuid2': a2, 'uuid3': a3})

        md_info = MetadataInfo(['uuid1', 'uuid2', 'uuid3'],
                               relative_fp='hmm.tsv')
        var = init_md_from_artifacts(md_info, ns, cfg)
        self.assertIsInstance(var, UsageVariable)
        self.assertEqual(var.var_type, 'metadata')
        rendered = var.use.render()
        exp = """from qiime2 import Metadata

thing1_a_0_md = thing1.view(Metadata)
thing2_a_0_md = thing2.view(Metadata)
thing3_a_0_md = thing3.view(Metadata)
merged_artifacts_0_md = thing1_a_0_md.merge(thing2_a_0_md, thing3_a_0_md)"""
        self.assertEqual(rendered, exp)

    def test_init_md_from_md_file_not_mdc(self):
        v0_uuid = '9f6a0f3e-22e6-4c39-8733-4e672919bbc7'
        v1_uuid = '0b8b47bd-f2f8-4029-923c-0e37a68340c3'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'):
            mixed = ProvDAG(os.path.join(DATA_DIR,
                            'mixed_v0_v1_uu_emperor.qzv'))
        v1_node = mixed.get_node_data(v1_uuid)
        md_id = 'whatevs'
        param_name = 'metadata'
        ns = UsageVarsDict({md_id: param_name})
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=pm)

        var = init_md_from_md_file(v1_node, param_name, md_id, ns, cfg)

        rendered = var.use.render()
        self.assertRegex(rendered, 'from qiime2 import Metadata')
        self.assertRegex(
            rendered,
            r'metadata_0_md = Metadata.load\(\<your metadata filepath\>\)')

    def test_init_md_from_md_file_md_is_mdc(self):
        dag = ProvDAG(os.path.join(DATA_DIR, 'v5_table.qza'))
        n_id = '99fa3670-aa1a-45f6-ba8e-803c976a1163'
        demux_node = dag.get_node_data(n_id)
        md_id = 'per_sample_sequences_0_barcodes'
        # We expect the variable name has already been added to the namespace
        ns = UsageVarsDict({md_id: 'barcodes'})
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=pm)

        var = init_md_from_md_file(demux_node, 'barcodes', md_id, ns, cfg)
        self.assertIsInstance(var, UsageVariable)
        self.assertEqual(var.var_type, 'column')
        rendered = var.use.render()
        mdc_name = 'barcodes_0_mdc_0'
        self.assertRegex(rendered, 'from qiime2 import Metadata')
        self.assertRegex(
            rendered,
            r"barcodes_0_md = Metadata.load\(<your metadata filepath>\)")
        self.assertRegex(rendered,
                         rf"{mdc_name} = barcodes_0_md."
                         r"get_column\(<column name>\)")


class BuildNoProvenanceUsageTests(CustomAssertions):
    def test_build_no_provenance_node_usage_w_complete_node(self):
        ns = NamespaceCollections()
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=pm)
        v0_uuid = '0b8b47bd-f2f8-4029-923c-0e37a68340c3'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'):
            single_no_prov = ProvDAG(os.path.join(DATA_DIR,
                                     'v0_uu_emperor.qzv'))
        v0 = single_no_prov.get_node_data(v0_uuid)
        build_no_provenance_node_usage(v0, v0_uuid, ns, cfg)
        out_var_name = '<visualization_0>'
        self.assertEqual(ns.usg_var_namespace, {v0_uuid: out_var_name})
        rendered = cfg.use.render()
        # Confirm the initial context comment is present once.
        self.assertREAppearsOnlyOnce(rendered, 'nodes have no provenance')
        header = '# Original Node ID                       String Description'
        self.assertREAppearsOnlyOnce(rendered, header)

        # Confirm expected values have been rendered
        exp_v0 = f'# {v0_uuid}   {out_var_name}'
        self.assertRegex(rendered, exp_v0)

    def test_build_no_provenance_node_usage_uuid_only_node(self):
        ns = NamespaceCollections()
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=pm)
        v0_uuid = '9f6a0f3e-22e6-4c39-8733-4e672919bbc7'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'):
            mixed = ProvDAG(os.path.join(DATA_DIR,
                            'mixed_v0_v1_uu_emperor.qzv'))
        # This node is only a parent UUID captured in the v1 node's action.yaml
        node = mixed.get_node_data(v0_uuid)
        self.assertEqual(node, None)
        build_no_provenance_node_usage(node, v0_uuid, ns, cfg)

        out_var_name = '<no-provenance-node_0>'
        self.assertEqual(ns.usg_var_namespace, {v0_uuid: out_var_name})

        rendered = cfg.use.render()
        # Confirm the initial context comment is present once.
        self.assertREAppearsOnlyOnce(rendered, 'nodes have no provenance')
        header = '# Original Node ID                       String Description'
        self.assertREAppearsOnlyOnce(rendered, header)

        # Confirm expected values have been rendered
        exp_v0 = f'# {v0_uuid}   {out_var_name.replace("-", "_")}'
        self.assertRegex(rendered, exp_v0)

    def test_build_no_provenance_node_usage_many_x(self):
        """
        Context should only be logged once.
        """
        ns = NamespaceCollections()
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=pm)

        # This function doesn't actually know about the DAG, so no need to join
        v0_uuid = '0b8b47bd-f2f8-4029-923c-0e37a68340c3'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'):
            single_no_prov = ProvDAG(os.path.join(DATA_DIR,
                                     'v0_uu_emperor.qzv'))
        v0_viz = single_no_prov.get_node_data(v0_uuid)

        dummy_viz_uuid = v0_uuid + '-dummy'
        # Will return the same type as v0, so will catch failure to deduplicate
        dummy_viz = single_no_prov.get_node_data(v0_uuid)

        tbl_uuid = '89af91c0-033d-4e30-8ac4-f29a3b407dc1'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{tbl_uuid}.*prior.*incomplete'):
            v0_tbl = ProvDAG(os.path.join(DATA_DIR, 'v0_table.qza'))
        tbl = v0_tbl.get_node_data(tbl_uuid)
        build_no_provenance_node_usage(v0_viz, v0_uuid, ns, cfg)
        build_no_provenance_node_usage(tbl, tbl_uuid, ns, cfg)
        build_no_provenance_node_usage(dummy_viz, dummy_viz_uuid, ns, cfg)
        self.assertIn(v0_uuid, ns.usg_var_namespace)
        self.assertIn(tbl_uuid, ns.usg_var_namespace)
        self.assertIn(dummy_viz_uuid, ns.usg_var_namespace)
        self.assertEqual(ns.usg_var_namespace[v0_uuid], '<visualization_0>')
        self.assertEqual(ns.usg_var_namespace[tbl_uuid],
                         '<feature_table_frequency_0>')
        self.assertEqual(ns.usg_var_namespace[dummy_viz_uuid],
                         '<visualization_1>')
        rendered = cfg.use.render()
        # Confirm the initial context isn't repeated.
        self.assertREAppearsOnlyOnce(rendered, 'nodes have no provenance')
        header = '# Original Node ID                       String Description'
        self.assertREAppearsOnlyOnce(rendered, header)

        # Confirm expected values have been rendered
        exp_v0 = '# 0b8b47bd-f2f8-4029-923c-0e37a68340c3   <visualization_0>'
        exp_t = ('# 89af91c0-033d-4e30-8ac4-f29a3b407dc1   '
                 '<feature_table_frequency_0>')
        exp_v1 = ('# 0b8b47bd-f2f8-4029-923c-0e37a68340c3-dummy   '
                  '<visualization_1>')
        self.assertRegex(rendered, exp_v0)
        self.assertRegex(rendered, exp_t)
        self.assertRegex(rendered, exp_v1)


class BuildImportUsageTests(CustomAssertions):
    def test_build_import_usage_python(self):
        ns = NamespaceCollections()
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=pm)
        dag = ProvDAG(os.path.join(DATA_DIR, 'v5_table.qza'))
        n_id = 'a35830e1-4535-47c6-aa23-be295a57ee1c'
        imp_node = dag.get_node_data(n_id)
        c_to_s_type = camel_to_snake(imp_node.type)
        unq_var_nm = c_to_s_type + '_0'
        build_import_usage(imp_node, ns, cfg)
        rendered = cfg.use.render()
        vars = ns.usg_vars
        out_name = vars[n_id].to_interface_name()

        self.assertIsInstance(vars[n_id], UsageVariable)
        self.assertEqual(vars[n_id].var_type, 'artifact')
        self.assertEqual(vars[n_id].name, unq_var_nm)
        self.assertRegex(rendered, 'from qiime2 import Artifact')
        self.assertRegex(rendered, rf'{out_name} = Artifact.import_data\(')
        self.assertRegex(rendered, imp_node.type)
        self.assertRegex(rendered, '<your data here>')

    def test_build_import_usage_cli(self):
        ns = NamespaceCollections()
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['cli'](),
                           use_recorded_metadata=False, pm=pm)
        dag = ProvDAG(os.path.join(DATA_DIR, 'v5_table.qza'))
        n_id = 'a35830e1-4535-47c6-aa23-be295a57ee1c'
        imp_node = dag.get_node_data(n_id)
        c_to_s_type = camel_to_snake(imp_node.type)
        unq_var_nm = c_to_s_type + '_0'
        build_import_usage(imp_node, ns, cfg)
        rendered = cfg.use.render()
        vars = ns.usg_vars
        out_name = vars[n_id].to_interface_name()

        self.assertIsInstance(vars[n_id], UsageVariable)
        self.assertEqual(vars[n_id].var_type, 'artifact')
        self.assertEqual(vars[n_id].name, unq_var_nm)
        self.assertRegex(rendered, r'qiime tools import \\')
        self.assertRegex(rendered, f"  --type '{imp_node.type}'")
        self.assertRegex(rendered, "  --input-path <your data here>")
        self.assertRegex(rendered, f"  --output-path {out_name}")


class BuildActionUsageTests(CustomAssertions):
    def test_build_action_usage_cli(self):
        plugin = 'demux'
        action = 'emp-single'
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['cli'](),
                           use_recorded_metadata=False, pm=pm)
        ns = NamespaceCollections()
        import_var = CLIUsageVariable(
            'imported_seqs_0', lambda: None, 'artifact', cfg.use)
        ns.usg_vars = {'a35830e1-4535-47c6-aa23-be295a57ee1c': import_var}
        dag = ProvDAG(os.path.join(DATA_DIR, 'v5_table.qza'))
        act_id = '3d69c8d1-a1fa-4ab3-ac88-3a98da15b2d5'
        n_id = '99fa3670-aa1a-45f6-ba8e-803c976a1163'
        node = dag.get_node_data(n_id)
        acts = ActionCollections(std_actions={act_id:
                                              {n_id: 'per_sample_sequences'}})
        unq_var_nm = node.action.output_name + '_0'
        build_action_usage(node, ns, acts.std_actions, act_id, cfg)
        rendered = cfg.use.render()
        out_name = ns.usg_vars[n_id].to_interface_name()

        vars = ns.usg_vars
        self.assertIsInstance(vars[n_id], UsageVariable)
        self.assertEqual(vars[n_id].var_type, 'artifact')
        self.assertEqual(vars[n_id].name, unq_var_nm)
        self.assertREAppearsOnlyOnce(rendered, "Replay attempts.*metadata")
        self.assertREAppearsOnlyOnce(rendered, "command may have received")
        act_undersc = re.sub('-', '_', action)
        self.assertREAppearsOnlyOnce(
            rendered,
            fr"saved at '.\/recorded_metadata\/{plugin}_{act_undersc}_0\/'")
        self.assertRegex(rendered, f"qiime {plugin} {action}")
        self.assertRegex(rendered, "--i-seqs imported-seqs-0.qza")
        self.assertRegex(rendered,
                         "--m-barcodes-file <your metadata filepath>")
        self.assertRegex(rendered, r"--m-barcodes-column <column name>")
        self.assertRegex(rendered, "--p-no-rev-comp-barcodes")
        self.assertRegex(rendered, "--p-no-rev-comp-mapping-barcodes")
        self.assertRegex(rendered, f"--o-per-sample-sequences {out_name}")

    def test_build_action_usage_cli_parameter_name_has_changed(self):
        plugin = 'emperor'
        action = 'plot'
        act_id = 'c147dfbc-139a-4db0-ac17-b11948247f93'
        pcoa_id = '9f6a0f3e-22e6-4c39-8733-4e672919bbc7'
        n_id = '0b8b47bd-f2f8-4029-923c-0e37a68340c3'
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['cli'](),
                           use_recorded_metadata=False, pm=pm)
        ns = NamespaceCollections()
        import_var = CLIUsageVariable(
            'pcoa', lambda: None, 'artifact', cfg.use)
        ns.usg_vars = {pcoa_id: import_var}
        mixed_uuid = '9f6a0f3e-22e6-4c39-8733-4e672919bbc7'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{mixed_uuid}.*prior.*incomplete'):
            dag = ProvDAG(os.path.join(DATA_DIR, 'mixed_v0_v1_uu_emperor.qzv'))
        node = dag.get_node_data(n_id)
        # This is a v1 node, so we don't have an output name. use type.
        out_name_raw = node.type.lower()
        acts = ActionCollections(std_actions={act_id:
                                              {n_id: out_name_raw}})
        unq_var_nm = out_name_raw + '_0'
        build_action_usage(node, ns, acts.std_actions, act_id, cfg)
        rendered = cfg.use.render()
        vars = ns.usg_vars
        out_name = vars[n_id].to_interface_name()

        self.assertIsInstance(vars[n_id], UsageVariable)
        self.assertEqual(vars[n_id].var_type, 'visualization')
        self.assertEqual(vars[n_id].name, unq_var_nm)

        self.assertRegex(rendered, f"qiime {plugin} {action}")
        self.assertRegex(rendered, "--i-pcoa pcoa.qza")
        self.assertRegex(rendered,
                         "--m-metadata-file <your metadata filepath>")
        self.assertRegex(rendered,
                         "(?s)parameter name was not found in your.*env")
        # This has become "custom-axes" since the .qzv was first recorded
        self.assertRegex(rendered, r"--\?-custom-axis DaysSinceExperimentSta")
        self.assertRegex(rendered, f"--o-visualization {out_name}")

    def test_build_action_usage_python(self):
        plugin = 'demux'
        action = 'emp_single'
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=pm)
        import_var = ArtifactAPIUsageVariable(
            'imported_seqs_0', lambda: None, 'artifact', cfg.use)
        seqs_id = 'a35830e1-4535-47c6-aa23-be295a57ee1c'
        ns = NamespaceCollections()
        ns.usg_vars = {seqs_id: import_var}
        dag = ProvDAG(os.path.join(DATA_DIR, 'v5_table.qza'))
        act_id = '3d69c8d1-a1fa-4ab3-ac88-3a98da15b2d5'
        n_id = '99fa3670-aa1a-45f6-ba8e-803c976a1163'
        node = dag.get_node_data(n_id)
        out_name_raw = node.action.output_name
        acts = ActionCollections(std_actions={act_id:
                                              {n_id: out_name_raw}})
        unq_var_nm = out_name_raw + '_0'
        build_action_usage(node, ns, acts.std_actions, act_id, cfg)
        rendered = cfg.use.render()
        vars = ns.usg_vars
        out_name = vars[n_id].to_interface_name()

        self.assertIsInstance(vars[n_id], UsageVariable)
        self.assertEqual(vars[n_id].var_type, 'artifact')
        self.assertEqual(vars[n_id].name, unq_var_nm)
        self.assertRegex(rendered, "from qiime2 import Metadata")
        self.assertRegex(
            rendered, f"import.*{plugin}.actions as {plugin}_actions")
        self.assertREAppearsOnlyOnce(rendered, "Replay attempts.*metadata")
        self.assertREAppearsOnlyOnce(rendered, "command may have received")
        self.assertREAppearsOnlyOnce(
            rendered,
            fr"saved at '.\/recorded_metadata\/{plugin}_{action}_0\/'")
        self.assertREAppearsOnlyOnce(rendered, "NOTE:.*substitute.*Metadata")
        md_name = 'barcodes_0_md'
        mdc_name = 'barcodes_0_mdc_0'
        self.assertRegex(rendered, rf'{md_name} = Metadata.load\(<.*filepath>')
        self.assertRegex(rendered,
                         rf'{mdc_name} = {md_name}.get_column\(<col')
        self.assertRegex(rendered,
                         rf'{out_name}, _ = {plugin}_actions.{action}\(')
        self.assertRegex(rendered, f'seqs.*{vars[seqs_id].name}')
        self.assertRegex(rendered, f'barcodes.*{mdc_name}')
        self.assertRegex(rendered, 'rev_comp_barcodes.*False')
        self.assertRegex(rendered, 'rev_comp_mapping_barcodes.*False')

    def test_build_action_usage_recorded_md(self):
        plugin = 'emperor'
        action = 'plot'
        md_param = 'metadata'
        act_id = 'c147dfbc-139a-4db0-ac17-b11948247f93'
        pcoa_id = '9f6a0f3e-22e6-4c39-8733-4e672919bbc7'
        n_id = '0b8b47bd-f2f8-4029-923c-0e37a68340c3'
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=True, pm=pm)
        import_var = ArtifactAPIUsageVariable(
            'pcoa', lambda: None, 'artifact', cfg.use)
        ns = NamespaceCollections()
        ns.usg_vars = {pcoa_id: import_var}
        mixed_uuid = '9f6a0f3e-22e6-4c39-8733-4e672919bbc7'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{mixed_uuid}.*prior.*incomplete'):
            dag = ProvDAG(os.path.join(DATA_DIR, 'mixed_v0_v1_uu_emperor.qzv'))
        node = dag.get_node_data(n_id)
        # This is a v1 node, so we don't have an output name. use type.
        out_name_raw = node.type.lower()
        acts = ActionCollections(std_actions={act_id:
                                              {n_id: out_name_raw}})
        unq_var_nm = out_name_raw + '_0'
        build_action_usage(node, ns, acts.std_actions, act_id, cfg)
        rendered = cfg.use.render()
        vars = ns.usg_vars
        out_name = vars[n_id].to_interface_name()

        self.assertIsInstance(vars[n_id], UsageVariable)
        self.assertEqual(vars[n_id].var_type, 'visualization')
        self.assertEqual(vars[n_id].name, unq_var_nm)

        self.assertRegex(rendered, "from qiime2 import Metadata")
        self.assertRegex(
            rendered, f"import.*{plugin}.actions as {plugin}_actions")

        md_name = f'{md_param}_0_md'
        self.assertRegex(rendered,
                         rf'{md_name} = Metadata.load\(.*metadata_0.tsv')

        self.assertRegex(rendered,
                         rf'{out_name}, = {plugin}_actions.{action}\(')
        self.assertRegex(rendered, f'pcoa.*{vars[pcoa_id].name}')
        self.assertRegex(rendered, f'metadata.*{md_name}')
        self.assertRegex(rendered, "custom_axis='DaysSinceExperimentStart'")

    @patch('qiime2.core.archive.provenance_lib.replay.init_md_from_artifacts')
    def test_build_action_usage_md_from_artifacts(self, patch):
        act_id = '59798956-9261-40f3-b70f-fc5059e379f5'
        n_id = '75f035ac-33fb-4d1c-bdcd-63ae1d564056'
        sd_act_id = '7862c023-13d4-4cd5-a9db-c92507055d25'
        sd_id = 'a42ea02f-8c40-432c-9b88-e602f6cd3787'
        dag = ProvDAG(os.path.join(DATA_DIR, 'md_tabulated_from_art.qzv'))
        cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS['python3'](),
                           use_recorded_metadata=False, pm=pm)
        sd_var = ArtifactAPIUsageVariable(
            'sample_data_alpha_diversity_0', lambda: None, 'artifact', cfg.use)
        ns = NamespaceCollections()
        ns.usg_vars = {sd_id: sd_var}
        node = dag.get_node_data(n_id)
        out_name_raw = node.action.output_name
        acts = ActionCollections(std_actions={act_id: {n_id: out_name_raw},
                                              sd_act_id: {sd_id: 'smpl_data'}})
        build_action_usage(node, ns, acts.std_actions, act_id, cfg)
        patch.assert_called_once_with(
            MetadataInfo(
                input_artifact_uuids=['a42ea02f-8c40-432c-9b88-e602f6cd3787'],
                relative_fp='input.tsv'), ns, cfg)


class BibContentTests(unittest.TestCase):
    def test_contents(self):
        series_21 = {
            'year': ' 2010 ',
            'title': ' Data Structures for Statistical Computing in Python ',
            'pages': ' 51 -- 56 ',
            'editor': ' Stéfan van der Walt and Jarrod Millman ',
            'booktitle': ' Proceedings of the 9th Python in Science Conferen',
            'author': ' Wes McKinney ',
            'ENTRYTYPE': 'inproceedings',
            'ID': 'view|types:2021.2.0|pandas.core.series:Series|0'}

        df_20 = {
            'year': ' 2010 ',
            'title': ' Data Structures for Statistical Computing in Python ',
            'pages': ' 51 -- 56 ',
            'editor': ' Stéfan van der Walt and Jarrod Millman ',
            'booktitle': ' Proceedings of the 9th Python in Science Conferen',
            'author': ' Wes McKinney ',
            'ENTRYTYPE': 'inproceedings',
            'ID': 'view|types:2020.2.0|pandas.core.frame:DataFrame|0'}

        self.assertEqual(BibContent(series_21), BibContent(df_20))
        self.assertEqual(hash(BibContent(series_21)), hash(BibContent(df_20)))
        # Set membership because these objects are equal and hash-equal
        self.assertIn(BibContent(series_21), {BibContent(df_20)})


class CitationsTests(unittest.TestCase):
    def test_dedupe_citations(self):
        fn = os.path.join(DATA_DIR, 'dupes.bib')
        with open(fn) as bibtex_file:
            bib_db = bp.load(bibtex_file)
        deduped = dedupe_citations(bib_db.entries)
        # Dedupe by DOI will preserve only one of the biom.table entries
        # Dedupe by contents should preserve only one of the pandas entries
        self.assertEqual(len(deduped), 3)
        # Confirm each paper is present. The len assertion ensures one-to-one
        lower_keys = [entry['ID'].lower() for entry in deduped]
        self.assertTrue(any('framework' in key for key in lower_keys))
        self.assertTrue(any('biom' in key for key in lower_keys))
        self.assertTrue(any('pandas' in key for key in lower_keys))

    def test_dedupe_pandas(self):
        """
        No match on ID, framework, or DOI, but matching content should
        deduplicate these
        """
        series_21 = {
            'year': ' 2010 ',
            'title': ' Data Structures for Statistical Computing in Python ',
            'pages': ' 51 -- 56 ',
            'editor': ' Stéfan van der Walt and Jarrod Millman ',
            'booktitle': ' Proceedings of the 9th Python in Science Conferen',
            'author': ' Wes McKinney ',
            'ENTRYTYPE': 'inproceedings',
            'ID': 'view|types:2021.2.0|pandas.core.series:Series|0'}

        df_20 = {
            'year': ' 2010 ',
            'title': ' Data Structures for Statistical Computing in Python ',
            'pages': ' 51 -- 56 ',
            'editor': ' Stéfan van der Walt and Jarrod Millman ',
            'booktitle': ' Proceedings of the 9th Python in Science Conferen',
            'author': ' Wes McKinney ',
            'ENTRYTYPE': 'inproceedings',
            'ID': 'view|types:2020.2.0|pandas.core.frame:DataFrame|0'}

        deduped = dedupe_citations([series_21, df_20])
        self.assertEqual(len(deduped), 1)

    def test_dedupe_silva(self):
        """
        These similar publications should not be deduped by content filter
        """
        s0 = {
            'year': '2007',
            'volume': '35',
            'title': 'SILVA: a comprehensive online resource for quality '
                     'checked and aligned ribosomal RNA sequence data '
                     'compatible with ARB',
            'pages': '7188-7196',
            'number': '21',
            'journal': 'Nucleic Acids Res',
            'author': 'Pruesse, Elmar and Quast, Christian and Knittel, Katrin'
                      ' and Fuchs, Bernhard M and Ludwig, Wolfgang and Peplies'
                      ', Jorg and Glockner, Frank Oliver',
            'ENTRYTYPE': 'article',
            'ID': 'action|rescript:2020.6.0+3.g772294c|'
                  'method:parse_silva_taxonomy|0'}
        s1 = {
            'year': '2013',
            'volume': '41',
            'title': 'The SILVA ribosomal RNA gene database project: '
                     'improved data processing and web-based tools',
            'publisher': 'Oxford University Press',
            'pages': 'D590-6',
            'number': 'Database issue',
            'journal': 'Nucleic Acids Res',
            'author': 'Quast, Christian and Pruesse, Elmar and Yilmaz, Pelin '
                      'and Gerken, Jan and Schweer, Timmy and Yarza, Pablo and'
                      ' Peplies, Jorg and Glockner, Frank Oliver',
            'ENTRYTYPE': 'article',
            'ID': 'action|rescript:2020.6.0+3.g772294c|'
                  'method:parse_silva_taxonomy|1'}
        deduped = dedupe_citations([s0, s1])
        self.assertEqual(len(deduped), 2)

    def test_collect_citations_no_deduped(self):
        dag = ProvDAG(TEST_DATA['5']['qzv_fp'])
        exp_keys = {'framework|qiime2:2018.11.0|0',
                    'action|feature-table:2018.11.0|method:rarefy|0',
                    'view|types:2018.11.0|BIOMV210DirFmt|0',
                    'view|types:2018.11.0|biom.table:Table|0',
                    'plugin|dada2:2018.11.0|0',
                    'action|alignment:2018.11.0|method:mafft|0',
                    'action|diversity:2018.11.0|method:beta_phylogenetic|0',
                    'action|diversity:2018.11.0|method:beta_phylogenetic|1',
                    'action|diversity:2018.11.0|method:beta_phylogenetic|2',
                    'action|diversity:2018.11.0|method:beta_phylogenetic|3',
                    'action|diversity:2018.11.0|method:beta_phylogenetic|4',
                    'view|types:2018.11.0|BIOMV210Format|0',
                    'plugin|emperor:2018.11.0|0',
                    'plugin|emperor:2018.11.0|1',
                    'action|phylogeny:2018.11.0|method:fasttree|0',
                    'action|alignment:2018.11.0|method:mask|0',
                    }
        citations = collect_citations(dag, deduplicate=False)
        keys = set(citations.entries_dict.keys())
        self.assertEqual(len(keys), len(exp_keys))
        self.assertEqual(keys, exp_keys)

    def test_collect_deduped(self):
        v5_tbl = ProvDAG(os.path.join(DATA_DIR, 'v5_table.qza'))
        std_keys = {'framework|qiime2:2018.11.0|0',
                    'view|types:2018.11.0|BIOMV210DirFmt|0',
                    'view|types:2018.11.0|biom.table:Table|0',
                    'plugin|dada2:2018.11.0|0'}
        citations = collect_citations(v5_tbl, deduplicate=False)
        keys = set(citations.entries_dict.keys())
        self.assertEqual(len(keys), len(std_keys))
        self.assertEqual(keys, std_keys)

        citations = collect_citations(v5_tbl, deduplicate=True)
        keys = set(citations.entries_dict.keys())
        # Dedupe by DOI will drop one of the biom.table entries
        self.assertEqual(len(keys), 3)
        # We want to confirm each paper is present - it doesn't matter which
        # biom entry is dropped.
        lower_keys = [key.lower() for key in keys]
        self.assertTrue(any('framework' in key for key in lower_keys))
        self.assertTrue(any('dada2' in key for key in lower_keys))
        self.assertTrue(any('biom' in key for key in lower_keys))

    def test_collect_citations_no_prov(self):
        v0_uuid = '9f6a0f3e-22e6-4c39-8733-4e672919bbc7'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'):
            mixed = ProvDAG(os.path.join(DATA_DIR,
                            'mixed_v0_v1_uu_emperor.qzv'))
        exp_keys = set()
        citations = collect_citations(mixed)
        keys = set(citations.entries_dict.keys())
        self.assertEqual(len(keys), 0)
        self.assertEqual(keys, exp_keys)

    def test_replay_citations(self):
        dag = ProvDAG(TEST_DATA['5']['qzv_fp'])
        exp_keys = ['framework|qiime2:2018.11.0|0',
                    'action|feature-table:2018.11.0|method:rarefy|0',
                    'view|types:2018.11.0|BIOMV210DirFmt|0',
                    'plugin|dada2:2018.11.0|0',
                    'action|alignment:2018.11.0|method:mafft|0',
                    'action|diversity:2018.11.0|method:beta_phylogenetic|0',
                    'action|diversity:2018.11.0|method:beta_phylogenetic|1',
                    'action|diversity:2018.11.0|method:beta_phylogenetic|2',
                    'action|diversity:2018.11.0|method:beta_phylogenetic|3',
                    'action|diversity:2018.11.0|method:beta_phylogenetic|4',
                    'plugin|emperor:2018.11.0|0',
                    'plugin|emperor:2018.11.0|1',
                    'action|phylogeny:2018.11.0|method:fasttree|0',
                    'action|alignment:2018.11.0|method:mask|0',
                    ]
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'citations.bib'
            out_fn = str(out_fp)
            replay_citations(dag, out_fn)
            self.assertTrue(out_fp.is_file())
            with open(out_fn, 'r') as fp:
                written = fp.read()
                for key in exp_keys:
                    self.assertIn(key, written)

    def test_replay_citations_no_prov(self):
        v0_uuid = '9f6a0f3e-22e6-4c39-8733-4e672919bbc7'
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'):
            mixed = ProvDAG(os.path.join(DATA_DIR,
                            'mixed_v0_v1_uu_emperor.qzv'))
        exp = "No citations were registered"
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'citations.bib'
            out_fn = str(out_fp)
            replay_citations(mixed, out_fn)
            self.assertTrue(out_fp.is_file())
            with open(out_fn, 'r') as fp:
                written = fp.read()
                self.assertIn(exp, written)


class WriteReproducibilitySupplementTests(CustomAssertions):
    def test_replay_supplement_from_fp(self):
        """
        Do we get expected zipfile contents when the supplement is
        created from a filepath
        """
        in_fp = TEST_DATA['5']['qzv_fp']
        in_fn = str(in_fp)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'supplement.zip'
            out_fn = str(out_fp)
            replay_supplement(
                payload=in_fn,
                out_fp=out_fn,
            )

            self.assertTrue(out_fp.is_file())
            self.assertTrue(zipfile.is_zipfile(out_fp))

            exp = {'python3_replay.py',
                   'cli_replay.sh',
                   'citations.bib',
                   'recorded_metadata/',
                   'recorded_metadata/demux_emp_single_0/barcodes_0.tsv',
                   }

            with zipfile.ZipFile(out_fp, 'r') as myzip:
                namelist_set = set(myzip.namelist())
                for item in exp:
                    self.assertIn(item, namelist_set)

    def test_replay_supplement_from_provdag(self):
        """
        Do we get expected zipfile contents when the supplement is
        created from a ProvDAG
        """
        in_fp = TEST_DATA['5']['qzv_fp']
        in_fn = str(in_fp)

        dag = ProvDAG(artifact_data=in_fn)

        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'supplement.zip'
            out_fn = str(out_fp)
            replay_supplement(
                payload=dag,
                out_fp=out_fn,
            )

            self.assertTrue(out_fp.is_file())
            self.assertTrue(zipfile.is_zipfile(out_fp))

            exp = {'python3_replay.py',
                   'cli_replay.sh',
                   'citations.bib',
                   'recorded_metadata/',
                   'recorded_metadata/demux_emp_single_0/barcodes_0.tsv',
                   }

            with zipfile.ZipFile(out_fp, 'r') as myzip:
                namelist_set = set(myzip.namelist())
                for item in exp:
                    self.assertIn(item, namelist_set)
