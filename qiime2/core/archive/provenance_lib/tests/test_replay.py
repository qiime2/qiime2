# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import bibtexparser as bp
import networkx as nx
import os
import pathlib
import shutil
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from qiime2 import Artifact
from qiime2.sdk import PluginManager
from qiime2.sdk.usage import Usage, UsageVariable
from qiime2.plugins import ArtifactAPIUsageVariable

from ..parse import ProvDAG
from ..replay import (
    ActionCollections, BibContent, ReplayConfig, ReplayNamespaces,
    UsageVariableRecord,
    build_no_provenance_node_usage, build_import_usage, build_action_usage,
    build_usage_examples, collect_citations, dedupe_citations,
    dump_recorded_md_file, group_by_action, init_md_from_artifacts,
    init_md_from_md_file, init_md_from_recorded_md, replay_provenance,
    uniquify_action_name, replay_citations
)
from .testing_utilities import CustomAssertions, DummyArtifacts
from ..usage_drivers import ReplayPythonUsage
from ...provenance import MetadataInfo

from qiime2.sdk.util import camel_to_snake


class ReplayNamespacesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.das = DummyArtifacts()
        cls.tempdir = cls.das.tempdir

    @classmethod
    def tearDownClass(cls):
        cls.das.free()

    # TODO
    def test_make_unique_name(self):
        pass

    # TODO
    def test_get_usg_var_uuid(self):
        pass

    def test_add_usage_var_workflow(self):
        """
        Smoke tests a common workflow with this data structure
        - Create a unique variable name by adding to .usg_var_namespace
        - Create a UsageVariable with that name
        - use the name to get the UUID (when we have Results, we have no UUIDs)
        - add the correctly-named UsageVariable to .usg_vars
        """
        use = Usage()
        uuid = self.das.concated_ints.uuid
        base_name = 'concated_ints'
        exp_name = base_name + '_0'
        ns = ReplayNamespaces()
        ns.add_usg_var_record(uuid, base_name)
        self.assertEqual(ns.get_usg_var_record(uuid).name, exp_name)

        def factory():  # pragma: no cover
            return Artifact.load(self.das.concated_ints.filepath)

        u_var = use.init_artifact(ns.get_usg_var_record(uuid).name, factory)
        self.assertEqual(u_var.name, exp_name)

        actual_uuid = ns.get_usg_var_uuid(u_var.name)
        self.assertEqual(actual_uuid, uuid)

        ns.update_usg_var_record(uuid, u_var)
        self.assertIsInstance(
            ns.get_usg_var_record(uuid).variable, UsageVariable
        )
        self.assertEqual(ns.get_usg_var_record(uuid).name, exp_name)


class ReplayProvenanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.das = DummyArtifacts()
        cls.tempdir = cls.das.tempdir

    @classmethod
    def tearDownClass(cls):
        cls.das.free()

    def test_replay_from_fp(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'rendered.txt'
            out_fn = str(out_fp)
            in_fn = self.das.concated_ints_with_md.filepath
            replay_provenance(
                ReplayPythonUsage, in_fn, out_fn, md_out_dir=tmpdir
            )

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
        in_fp = self.das.concated_ints.filepath
        with self.assertRaisesRegex(
                ValueError, "Metadata not parsed for replay. Re-run"
        ):
            replay_provenance(
                ReplayPythonUsage,
                in_fp,
                'unused_fp',
                parse_metadata=False,
                use_recorded_metadata=True
            )

    def test_replay_dump_md_without_parse(self):
        in_fp = self.das.concated_ints.filepath
        with self.assertRaisesRegex(
                ValueError, "(?s)Metadata not parsed,.*dump_recorded_meta"
        ):
            replay_provenance(
                ReplayPythonUsage,
                in_fp,
                'unused_fp',
                parse_metadata=False,
                dump_recorded_metadata=True
            )

    def test_replay_md_out_dir_without_parse(self):
        in_fp = self.das.concated_ints.filepath
        with self.assertRaisesRegex(
                ValueError, "(?s)Metadata not parsed,.*not.*metadata output"
        ):
            replay_provenance(
                ReplayPythonUsage,
                in_fp,
                'unused_fp',
                parse_metadata=False,
                dump_recorded_metadata=False,
                md_out_dir='/user/dumb/some_filepath'
            )

    def test_replay_use_md_without_dump_md(self):
        in_fp = self.das.concated_ints.filepath
        with self.assertRaisesRegex(
                NotImplementedError,
                "(?s)uses.*metadata.*must.*written to disk"
        ):
            replay_provenance(
                ReplayPythonUsage,
                in_fp,
                'unused_fp',
                use_recorded_metadata=True,
                dump_recorded_metadata=False
            )

    def test_replay_from_provdag(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'rendered.txt'
            out_fn = str(out_fp)
            dag = self.das.concated_ints_with_md.dag
            replay_provenance(ReplayPythonUsage, dag, out_fn,
                              md_out_dir=tmpdir)

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
        dag = ProvDAG(self.das.concated_ints_with_md.filepath,
                      validate_checksums=False,
                      parse_metadata=False)
        with self.assertRaisesRegex(
                ValueError, "Metadata not parsed for replay"
        ):
            replay_provenance(
                ReplayPythonUsage, dag, 'unused', use_recorded_metadata=True
            )

    def test_replay_from_provdag_ns_collision(self):
        """
        This artifact's dag contains a few results with the output-name
        filtered-table, so is a good check for namespace collisions if
        we're not uniquifying variable names properly.
        """
        with tempfile.TemporaryDirectory() as tempdir:
            self.das.concated_ints.artifact.save(
                os.path.join(tempdir, 'c1.qza'))
            self.das.other_concated_ints.artifact.save(
                os.path.join(tempdir, 'c2.qza'))
            dag = ProvDAG(tempdir)

        exp = ['concatenated_ints_0', 'concatenated_ints_1']
        with tempfile.TemporaryDirectory() as tempdir:
            out_path = pathlib.Path(tempdir) / 'ns_coll.txt'
            replay_provenance(
                 ReplayPythonUsage, dag, out_path, md_out_dir=tempdir
            )

            with open(out_path, 'r') as fp:
                rendered = fp.read()
                for name in exp:
                    self.assertIn(name, rendered)

    def test_replay_optional_param_is_none(self):
        dag = self.das.int_seq_optional_input.dag
        with tempfile.TemporaryDirectory() as tempdir:
            out_path = pathlib.Path(tempdir) / 'ns_coll.txt'

            replay_provenance(
                ReplayPythonUsage, dag, out_path, md_out_dir=tempdir
            )
            with open(out_path, 'r') as fp:
                rendered = fp.read()
            self.assertIn('ints=int_sequence1_0', rendered)
            self.assertIn('num1=', rendered)
            self.assertNotIn('optional1=', rendered)
            self.assertNotIn('num2=', rendered)


class MultiplePluginTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from qiime2.sdk.plugin_manager import PluginManager
        from qiime2 import Artifact

        cls.pm = PluginManager()
        cls.dp = cls.pm.plugins['dummy-plugin']
        cls.op = cls.pm.plugins['other-plugin']
        cls.tempdir = tempfile.mkdtemp(prefix='qiime2-other-plugin-temp-')

        int_seq = Artifact.import_data('IntSequence1', [1, 2, 3, 4])
        concat_ints = cls.op.methods['concatenate_ints']
        split_ints = cls.dp.methods['split_ints']

        concated_ints, = concat_ints(
            int_seq, int_seq, int_seq, 5, 6
        )
        cls.splitted_ints, _ = split_ints(concated_ints)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tempdir)

    def test_multiple_plugins_in_provenance(self):
        fp = os.path.join(self.tempdir, 'splitted_ints.qza')
        self.splitted_ints.save(fp)

        with tempfile.TemporaryDirectory() as tempdir:
            out_fp = os.path.join(tempdir, 'rendered.txt')
            replay_provenance(
                ReplayPythonUsage, fp, out_fp, md_out_dir=tempdir
            )

            with open(out_fp, 'r') as fp:
                rendered = fp.read()

        self.assertIn('from qiime2 import Artifact', rendered)
        self.assertIn(
            'import qiime2.plugins.dummy_plugin.actions as '
            'dummy_plugin_actions',
            rendered
        )
        self.assertIn(
            'import qiime2.plugins.other_plugin.actions as '
            'other_plugin_actions',
            rendered
        )
        self.assertIn(
            'dummy_plugin_actions.split_ints(', rendered
        )
        self.assertIn(
            'other_plugin_actions.concatenate_ints(', rendered
        )


class ReplayProvDAGDirectoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.das = DummyArtifacts()
        cls.tempdir = cls.das.tempdir

    @classmethod
    def tearDownClass(cls):
        cls.das.free()

    def test_directory_replay_multiple_imports(self):
        """
        The directory being parsed here contains two pairs of duplicates,
        and should replay as only two import statements.
        """
        outer_dir = os.path.join(self.tempdir, 'outer')
        inner_dir = os.path.join(outer_dir, 'inner')
        os.makedirs(inner_dir)

        for artifact in self.das.single_int, self.das.single_int2:
            for dir_ in inner_dir, outer_dir:
                artifact.artifact.save(
                    os.path.join(dir_, f'{artifact.name}.qza')
                )

        dir_dag = ProvDAG(outer_dir)
        self.assertEqual(len(dir_dag._parsed_artifact_uuids), 2)
        self.assertIn(self.das.single_int.uuid, dir_dag.dag)
        self.assertIn(self.das.single_int2.uuid, dir_dag.dag)

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
            replay_provenance(
                ReplayPythonUsage, dir_dag, out_path, md_out_dir=tempdir
            )
            self.assertTrue(out_path.is_file())

            with open(out_path, 'r') as fp:
                rendered = fp.read()
                self.assertRegex(rendered, exp_1)
                self.assertRegex(rendered, exp_2)


class BuildUsageExamplesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.das = DummyArtifacts()
        cls.tempdir = cls.das.tempdir
        cls.pm = PluginManager()

    @classmethod
    def tearDownClass(cls):
        cls.das.free()

    @patch('qiime2.core.archive.provenance_lib.replay.build_action_usage')
    @patch('qiime2.core.archive.provenance_lib.replay.build_import_usage')
    @patch('qiime2.core.archive.provenance_lib.replay.'
           'build_no_provenance_node_usage')
    def test_build_usage_examples(self, n_p_builder, imp_builder, act_builder):
        ns = ReplayNamespaces()
        dag = self.das.concated_ints_with_md.dag
        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False,
            pm=self.pm
        )
        build_usage_examples(dag, cfg, ns)

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
        ns = ReplayNamespaces()
        uuid = self.das.table_v0.uuid
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{uuid}.*prior.*incomplete'
        ):
            dag = ProvDAG(self.das.table_v0.filepath)

        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False,
            pm=self.pm
        )
        build_usage_examples(dag, cfg, ns)

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
        shutil.copy(self.das.table_v0.filepath, mixed_dir)
        shutil.copy(self.das.concated_ints_v6.filepath, mixed_dir)

        ns = ReplayNamespaces()
        v0_uuid = self.das.table_v0.uuid
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'
        ):
            dag = ProvDAG(mixed_dir)

        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False,
            pm=self.pm
        )
        build_usage_examples(dag, cfg, ns)

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
        shutil.copy(self.das.concated_ints_with_md.filepath, many_dir)
        shutil.copy(self.das.splitted_ints.filepath, many_dir)
        shutil.copy(self.das.pipeline_viz.filepath, many_dir)

        ns = ReplayNamespaces()
        dag = ProvDAG(many_dir)
        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False,
            pm=self.pm
        )
        build_usage_examples(dag, cfg, ns)

        n_p_builder.assert_not_called()
        # concated_ints_with_md is loaded from disk so imports don't overlap
        # with splitted_ints and pipeline_viz
        self.assertEqual(imp_builder.call_count, 6)
        self.assertEqual(act_builder.call_count, 4)


class MiscHelperFnTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.das = DummyArtifacts()
        cls.tempdir = cls.das.tempdir
        cls.pm = PluginManager()

    @classmethod
    def tearDownClass(cls):
        cls.das.free()

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

    def test_dump_recorded_md_file_no_md(self):
        uuid = self.das.table_v0.uuid
        dag = self.das.table_v0.dag

        cfg = ReplayConfig(use=ReplayPythonUsage(),
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
    @classmethod
    def setUpClass(cls):
        cls.das = DummyArtifacts()
        cls.tempdir = cls.das.tempdir
        cls.pm = PluginManager()

    @classmethod
    def tearDownClass(cls):
        cls.das.free()

    def test_gba_with_provenance(self):
        self.maxDiff = None

        ns = ReplayNamespaces()
        dag = self.das.concated_ints_v6.dag
        sorted_nodes = nx.topological_sort(dag.collapsed_view)
        actual = group_by_action(dag, sorted_nodes, ns)
        exp = {
            'b49e497c-19b2-49f7-b9a2-0d837016c151': {
                '8dea2f1a-2164-4a85-9f7d-e0641b1db22b': 'int_sequence1'
            },
            '12988290-1ebf-47ad-8c34-5469d42e5ffe': {
                '7727c060-5384-445d-b007-b64b41a090ee': 'int_sequence2'
            },
            '5035a60e-6f9a-40d4-b412-48ae52255bb5': {
                '6facaf61-1676-45eb-ada0-d530be678b27': 'concatenated_ints'
            }
        }
        self.assertEqual(actual.std_actions, exp)
        self.assertEqual(actual.no_provenance_nodes, [])

    def test_gba_no_provenance(self):
        ns = ReplayNamespaces()
        dag = self.das.table_v0.dag
        uuid = self.das.table_v0.uuid

        sorted_nodes = nx.topological_sort(dag.collapsed_view)
        action_collections = group_by_action(dag, sorted_nodes, ns)
        self.assertEqual(action_collections.std_actions, {})
        self.assertEqual(action_collections.no_provenance_nodes, [uuid])

    def test_gba_some_nodes_missing_provenance(self):
        mixed_dir = os.path.join(self.tempdir, 'mixed-dir')
        os.mkdir(mixed_dir)
        shutil.copy(self.das.table_v0.filepath, mixed_dir)
        shutil.copy(self.das.concated_ints_v6.filepath, mixed_dir)

        ns = ReplayNamespaces()
        v0_uuid = self.das.table_v0.uuid
        with self.assertWarnsRegex(
                UserWarning, f'(:?)Art.*{v0_uuid}.*prior.*incomplete'
        ):
            dag = ProvDAG(mixed_dir)

        sorted_nodes = nx.topological_sort(dag.collapsed_view)
        action_collections = group_by_action(dag, sorted_nodes, ns)

        exp = {
            'b49e497c-19b2-49f7-b9a2-0d837016c151': {
                '8dea2f1a-2164-4a85-9f7d-e0641b1db22b': 'int_sequence1'
            },
            '12988290-1ebf-47ad-8c34-5469d42e5ffe': {
                '7727c060-5384-445d-b007-b64b41a090ee': 'int_sequence2'
            },
            '5035a60e-6f9a-40d4-b412-48ae52255bb5': {
                '6facaf61-1676-45eb-ada0-d530be678b27': 'concatenated_ints'
            }
        }
        self.assertEqual(action_collections.std_actions, exp)
        self.assertEqual(action_collections.no_provenance_nodes, [v0_uuid])


class InitializerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.das = DummyArtifacts()
        cls.tempdir = cls.das.tempdir
        cls.pm = PluginManager()

        with zipfile.ZipFile(cls.das.concated_ints_with_md.filepath) as zf:
            root_node_id = cls.das.concated_ints_with_md.uuid
            all_filenames = zf.namelist()
            dag = cls.das.concated_ints_with_md.dag
            for node in dag.nodes:
                md_path = os.path.join(
                    root_node_id, 'provenance', 'artifacts', node, 'action',
                    'metadata.tsv'
                )
                if md_path in all_filenames:
                    cls.md_node_id = node
                else:
                    cls.non_md_node_id = node

        with zipfile.ZipFile(
            cls.das.concated_ints_with_md_column.filepath
        ) as zf:
            root_node_id = cls.das.concated_ints_with_md_column.uuid
            all_filenames = zf.namelist()
            dag = cls.das.concated_ints_with_md_column.dag
            for node in dag.nodes:
                md_path = os.path.join(
                    root_node_id, 'provenance', 'artifacts', node, 'action',
                    'metadata.tsv'
                )
                if md_path in all_filenames:
                    cls.mdc_node_id = node
                else:
                    cls.non_mdc_node_id = node

    @classmethod
    def tearDownClass(cls):
        cls.das.free()

    def test_init_md_from_artifacts_no_artifacts(self):
        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False,
            pm=self.pm
        )
        ns = ReplayNamespaces

        # create dummy hash '0', not relevant here
        md_info = MetadataInfo([], 'hmm.tsv', '0')
        with self.assertRaisesRegex(
            ValueError, "not.*used.*input_artifact_uuids.*empty"
        ):
            init_md_from_artifacts(md_info, ns, cfg)

    def test_init_md_from_artifacts_one_art(self):
        # This helper doesn't capture real data, so we're only smoke testing,
        # checking type, and confirming the repr looks reasonable.
        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False,
            pm=self.pm
        )

        # We expect artifact vars have already been added to the namespace
        a1 = cfg.use.init_artifact(name='thing1', factory=lambda: None)
        ns = ReplayNamespaces()
        ns._usg_var_ns = {'uuid1': UsageVariableRecord('thing1', a1)}

        # create dummy hash '0', not relevant here
        md_info = MetadataInfo(['uuid1'], 'hmm.tsv', '0')
        var = init_md_from_artifacts(md_info, ns, cfg)
        self.assertIsInstance(var, UsageVariable)
        self.assertEqual(var.var_type, 'metadata')
        rendered = var.use.render()

        self.assertIn('from qiime2 import Metadata', rendered)
        self.assertIn('thing1_a_0_md = thing1.view(Metadata)', rendered)

    def test_init_md_from_artifacts_many(self):
        # This helper doesn't capture real data, so we're only smoke testing,
        # checking type, and confirming the repr looks reasonable.
        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False,
            pm=self.pm
        )

        # We expect artifact vars have already been added to the namespace
        a1 = cfg.use.init_artifact(name='thing1', factory=lambda: None)
        a2 = cfg.use.init_artifact(name='thing2', factory=lambda: None)
        a3 = cfg.use.init_artifact(name='thing3', factory=lambda: None)

        ns = ReplayNamespaces()
        ns._usg_var_ns = {
            'uuid1': UsageVariableRecord('thing1', a1),
            'uuid2': UsageVariableRecord('thing2', a2),
            'uuid3': UsageVariableRecord('thing3', a3),
        }

        # create dummy hash '0', not relevant here
        md_info = MetadataInfo(['uuid1', 'uuid2', 'uuid3'], 'hmm.tsv', '0')
        var = init_md_from_artifacts(md_info, ns, cfg)
        self.assertIsInstance(var, UsageVariable)
        self.assertEqual(var.var_type, 'metadata')
        rendered = var.use.render()

        self.assertIn('from qiime2 import Metadata', rendered)
        self.assertIn('thing1_a_0_md = thing1.view(Metadata)', rendered)
        self.assertIn('thing2_a_0_md = thing2.view(Metadata)', rendered)
        self.assertIn('thing3_a_0_md = thing3.view(Metadata)', rendered)
        self.assertIn(
            'merged_artifacts_0_md = '
            'thing1_a_0_md.merge(thing2_a_0_md, thing3_a_0_md)',
            rendered
        )

    def test_init_md_from_md_file(self):
        dag = self.das.concated_ints_with_md.dag
        md_node = dag.get_node_data(self.md_node_id)
        md_id = 'whatevs'
        param_name = 'metadata'

        ns = ReplayNamespaces()
        ns.add_usg_var_record(md_id, param_name)

        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False,
            pm=self.pm
        )

        var = init_md_from_md_file(md_node, param_name, md_id, ns, cfg)

        rendered = var.use.render()
        self.assertIn('from qiime2 import Metadata', rendered)
        self.assertIn(
            'metadata_0_md = Metadata.load(<your metadata filepath>)',
            rendered
        )

    def test_init_md_from_recorded_md(self):
        dag = self.das.concated_ints_with_md.dag
        no_md_node = dag.get_node_data(self.non_md_node_id)
        md_node = dag.get_node_data(self.md_node_id)
        var_name = 'metadata_0'
        param_name = 'metadata'

        ns = ReplayNamespaces()
        ns.add_usg_var_record(var_name, param_name)

        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False,
            pm=self.pm
        )

        md_fn = 'identity_with_metadata/metadata_0'

        with self.assertRaisesRegex(ValueError, 'only.*call.*if.*metadata'):
            init_md_from_recorded_md(
                no_md_node, param_name, var_name, ns, cfg, md_fn
            )

        var = init_md_from_recorded_md(
            md_node, param_name, var_name, ns, cfg, md_fn
        )
        self.assertIsInstance(var, UsageVariable)
        self.assertEqual(var.var_type, 'metadata')

        rendered = cfg.use.render()
        self.assertIn('from qiime2 import Metadata', rendered)
        self.assertIn('metadata_0_md = Metadata.load', rendered)
        self.assertIn(
            'recorded_metadata/identity_with_metadata/metadata_0', rendered
        )

    def test_init_md_from_recorded_mdc(self):
        dag = self.das.concated_ints_with_md_column.dag
        no_md_node = dag.get_node_data(self.non_mdc_node_id)
        md_node = dag.get_node_data(self.mdc_node_id)
        var_name = 'metadata_0'
        param_name = 'metadata'

        ns = ReplayNamespaces()
        ns.add_usg_var_record(var_name, param_name)

        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False,
            pm=self.pm
        )

        md_fn = 'identity_with_metadata_column/metadata_0'

        with self.assertRaisesRegex(ValueError, 'only.*call.*if.*metadata'):
            init_md_from_recorded_md(
                no_md_node, param_name, var_name, ns, cfg, md_fn
            )

        var = init_md_from_recorded_md(
            md_node, param_name, var_name, ns, cfg, md_fn
        )
        self.assertIsInstance(var, UsageVariable)
        self.assertEqual(var.var_type, 'column')

        rendered = cfg.use.render()
        self.assertIn('from qiime2 import Metadata', rendered)
        self.assertIn('metadata_0_md = Metadata.load', rendered)
        self.assertIn('.get_column(', rendered)
        self.assertIn('recorded_metadata/identity_with_metadata_column/'
                      'metadata_0.tsv', rendered)


class BuildNoProvenanceUsageTests(CustomAssertions):
    @classmethod
    def setUpClass(cls):
        cls.das = DummyArtifacts()
        cls.tempdir = cls.das.tempdir
        cls.pm = PluginManager()

    @classmethod
    def tearDownClass(cls):
        cls.das.free()

    def test_build_no_provenance_node_usage_w_complete_node(self):
        ns = ReplayNamespaces()
        cfg = ReplayConfig(use=ReplayPythonUsage(),
                           use_recorded_metadata=False, pm=self.pm)
        uuid = self.das.table_v0.uuid
        dag = self.das.table_v0.dag
        v0_node = dag.get_node_data(uuid)
        build_no_provenance_node_usage(v0_node, uuid, ns, cfg)

        out_var_name = 'feature_table_frequency_0'
        self.assertIn(uuid, ns._usg_var_ns)
        self.assertEqual(ns._usg_var_ns[uuid].name, out_var_name)

        rendered = cfg.use.render()
        # Confirm the initial context comment is present once.
        self.assertREAppearsOnlyOnce(rendered, 'nodes have no provenance')
        header = '# Original Node ID                       String Description'
        self.assertREAppearsOnlyOnce(rendered, header)

        # Confirm expected values have been rendered
        exp_v0 = f'# {uuid}   feature_table_frequency_0'
        self.assertRegex(rendered, exp_v0)

    def test_build_no_provenance_node_usage_uuid_only_node(self):
        ns = ReplayNamespaces()
        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False,
            pm=self.pm
        )

        uuid = 'some-uuid'
        node = None
        build_no_provenance_node_usage(node, uuid, ns, cfg)

        out_var_name = 'no-provenance-node_0'
        self.assertIn(uuid, ns._usg_var_ns)
        self.assertEqual(ns._usg_var_ns[uuid].name, out_var_name)

        rendered = cfg.use.render()
        # Confirm the initial context comment is present once.
        self.assertREAppearsOnlyOnce(rendered, 'nodes have no provenance')
        header = '# Original Node ID                       String Description'
        self.assertREAppearsOnlyOnce(rendered, header)

        # Confirm expected values have been rendered
        exp_v0 = f'# {uuid}   no_provenance_node_0'
        self.assertRegex(rendered, exp_v0)

    def test_build_no_provenance_node_usage_many(self):
        ns = ReplayNamespaces()
        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False, pm=self.pm
        )

        # This function doesn't actually know about the DAG, so no need to join
        uuid = self.das.table_v0.uuid
        dag = self.das.table_v0.dag
        v0_node = dag.get_node_data(uuid)

        dummy_node_uuid = uuid + '-dummy'
        dummy_node = dag.get_node_data(uuid)

        build_no_provenance_node_usage(v0_node, uuid, ns, cfg)
        build_no_provenance_node_usage(dummy_node, dummy_node_uuid, ns, cfg)
        self.assertIn(uuid, ns._usg_var_ns)
        self.assertIn(dummy_node_uuid, ns._usg_var_ns)
        self.assertEqual(
            ns._usg_var_ns[uuid].name, 'feature_table_frequency_0'
        )
        self.assertEqual(
            ns._usg_var_ns[dummy_node_uuid].name, 'feature_table_frequency_1'
        )

        rendered = cfg.use.render()
        # Confirm the initial context isn't repeated.
        self.assertREAppearsOnlyOnce(rendered, 'nodes have no provenance')
        header = '# Original Node ID                       String Description'
        self.assertREAppearsOnlyOnce(rendered, header)

        # Confirm expected values have been rendered
        exp_og = f'# {uuid}   feature_table_frequency_0'
        exp_dummy = f'# {uuid}-dummy   feature_table_frequency_1'
        self.assertRegex(rendered, exp_og)
        self.assertRegex(rendered, exp_dummy)


class BuildImportUsageTests(CustomAssertions):
    @classmethod
    def setUpClass(cls):
        cls.das = DummyArtifacts()
        cls.tempdir = cls.das.tempdir
        cls.pm = PluginManager()

    @classmethod
    def tearDownClass(cls):
        cls.das.free()

    def test_build_import_usage_python(self):
        ns = ReplayNamespaces()
        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False,
            pm=self.pm
        )

        dag = self.das.concated_ints_v6.dag
        import_uuid = '8dea2f1a-2164-4a85-9f7d-e0641b1db22b'
        import_node = dag.get_node_data(import_uuid)
        c_to_s_type = camel_to_snake(import_node.type)
        unq_var_nm = c_to_s_type + '_0'
        build_import_usage(import_node, ns, cfg)


        usg_var = ns.get_usg_var_record(import_uuid).variable
        self.assertIsInstance(usg_var, UsageVariable)
        self.assertEqual(usg_var.var_type, 'artifact')
        self.assertEqual(usg_var.name, unq_var_nm)

        rendered = cfg.use.render()
        out_name = usg_var.to_interface_name()
        self.assertRegex(rendered, 'from qiime2 import Artifact')
        self.assertRegex(rendered, rf'{out_name} = Artifact.import_data\(')
        self.assertRegex(rendered, import_node.type)
        self.assertRegex(rendered, '<your data here>')


class ReplayResultCollectionTests(CustomAssertions):
    '''
    One of three structures may be encoutered when parsing the inputs section
    of action.yaml, described below:

    case 1 (single artifact case):

        inputs:
        - some_input_name: some_uuid
        - some_other_input_name: some_other_uuid
        (...)

        For the single artifact case we may find the artifact in the usage
        variable namespace or it might only exist in a result collection.
        - case 1a:
            The artifact exists in the usage variable namespace so access it
            from there, not from any result collection that may contain it.
        - case 1b:
            The artifact does not exist in the usage variable namespace so
            find it in a result collection, destructure it, and then use the
            destructured artifact.

    case 2 (list of artifacts case):

        inputs:
        - some_input_name:
            - some_uuid
            - some_other_uuid
        (...)

        For the list of artifacts case the list contents may be equivalent to
        an existing result collection, or it may not be.
        - case 2a:
            The list contents are not equivalent to any existing result
            collection so for each member do either case 1a or 1b as described
            above.
        - case 2b:
            The list contents are equivalent to an existing result collection
            so pass in the result collection directly (which will be cast to
            a list).

    case 3 (result collection case):

        inputs:
        - result_collection_name:
            - some_key: some_uuid
            - some_other_key: some_other_uuid
        (...)

        For the result collection case an equivalent existing result collection
        may exist or it may not.
        - case 3a:
            No equivalent result collection is found so for each member follow
            case 1a or 1b as described above.
        - case 3b:
            An equivalent result collection is found so pass it in directly.

    '''
    @classmethod
    def setUpClass(cls):
        cls.pm = PluginManager()
        cls.dp = cls.pm.plugins['dummy-plugin']

        cls.single_int = Artifact.import_data('SingleInt', 0)
        cls.dict_of_ints = cls.dp.methods['dict_of_ints']
        cls.list_of_ints = cls.dp.methods['list_of_ints']

    def test_cases_1a_1b(self):
        '''
        The `single_int` usage variable is not found before the first call to
        `optional_artifact_pipeline`, so it is accessed from the result
        collection it belongs to (case 1b) and added to the usage variable
        namespace. When we call `optional_artifact_pipeline` again it should
        not destructure the result collection a second time but instead just
        use the available usage variable that we added to the namespace
        previously (case 1a).
        '''
        int_seq = Artifact.import_data('IntSequence1', [1, 1, 2])

        opt_art_pipeline = self.dp.pipelines['optional_artifact_pipeline']

        rc = {'int1': self.single_int}

        rc_out, = self.dict_of_ints(rc)
        int_destructured = rc_out['int1']
        int_seq_out, = opt_art_pipeline(int_seq, int_destructured)
        int_seq_out, = opt_art_pipeline(int_seq_out, int_destructured)

        with tempfile.TemporaryDirectory() as tempdir:
            in_fp = pathlib.Path(tempdir) / 'int_seq_out.qza'
            int_seq_out.save(in_fp)
            dag = ProvDAG(in_fp)

            out_fp = pathlib.Path(tempdir) / 'rendered.txt'
            out_fn = str(out_fp)
            replay_provenance(ReplayPythonUsage, dag, out_fn)

            with open(out_fp) as fh:
                rendered = fh.read()

        exp = '''\
single_int_1 = output_0_artifact_collection['int1']
ints_1, = dummy_plugin_actions.optional_artifact_pipeline(
    int_sequence=int_sequence1_0,
    single_int=single_int_1,
)
# SAVE: comment out the following with '# ' to skip saving Results to disk
ints_1.save('ints_1')

ints_2, = dummy_plugin_actions.optional_artifact_pipeline(
    int_sequence=ints_1,
    single_int=single_int_1,
)
'''
        self.assertIn(exp, rendered)

    def test_case_2a(self):
        rc = {'int1': self.single_int, 'int2': self.single_int}
        list_rc_out, = self.list_of_ints(rc)
        int1_destructured = list_rc_out['0']

        with tempfile.TemporaryDirectory() as tempdir:
            in_fp = pathlib.Path(tempdir) / 'int1_destructured.qza'
            int1_destructured.save(in_fp)
            dag = ProvDAG(in_fp)

            out_fp = pathlib.Path(tempdir) / 'rendered.txt'
            out_fn = str(out_fp)
            replay_provenance(ReplayPythonUsage, dag, out_fn)

            with open(out_fp) as fh:
                rendered = fh.read()

        exp = '''\
output_0_artifact_collection, = dummy_plugin_actions.list_of_ints(
    ints=[single_int_0, single_int_0],
)
'''
        self.assertIn(exp, rendered)

    def test_case_2b(self):
        rc = {'int1': self.single_int}
        rc_out, = self.dict_of_ints(rc)
        list_rc_out, = self.list_of_ints(rc_out)
        int1_destructured = list_rc_out['0']

        with tempfile.TemporaryDirectory() as tempdir:
            in_fp = pathlib.Path(tempdir) / 'int1_destructured.qza'
            int1_destructured.save(in_fp)
            dag = ProvDAG(in_fp)

            out_fp = pathlib.Path(tempdir) / 'rendered.txt'
            out_fn = str(out_fp)
            replay_provenance(ReplayPythonUsage, dag, out_fn)

            with open(out_fp) as fh:
                rendered = fh.read()

        exp = '''\
output_1_artifact_collection, = dummy_plugin_actions.list_of_ints(
    ints=output_0_artifact_collection,
)
'''
        self.assertIn(exp, rendered)

    def test_case_3a(self):
        '''
        The `int1_destructured` artifact won't be found in the usage variable
        namespace, while `self.single_int` will be, so assert that int1 is
        destructured and `self.single_int` is just used as is because its
        variable name is already available.
        '''
        rc1 = {'int1': self.single_int}
        rc1_out, = self.dict_of_ints(rc1)
        int1_destructured = rc1_out['int1']

        rc2 = {'int1': int1_destructured, 'int2': self.single_int}
        rc2_out, = self.dict_of_ints(rc2)
        int3_destructured = rc2_out['int2']

        with tempfile.TemporaryDirectory() as tempdir:
            in_fp = pathlib.Path(tempdir) / 'int2_destructured.qza'
            int3_destructured.save(in_fp)
            dag = ProvDAG(in_fp)

            out_fp = pathlib.Path(tempdir) / 'rendered.txt'
            out_fn = str(out_fp)
            replay_provenance(ReplayPythonUsage, dag, out_fn)

            with open(out_fp) as fh:
                rendered = fh.read()

        exp = '''\
ints_1 = output_0_artifact_collection['int1']
ints_2_artifact_collection = ResultCollection({
    'int1': ints_1,
    'int2': single_int_0,
})
output_1_artifact_collection, = dummy_plugin_actions.dict_of_ints(
    ints=ints_2_artifact_collection,
)
'''
        self.assertIn(exp, rendered)
        self.assertREAppearsOnlyOnce(rendered, r'\[.*\]')

    def test_case_3b(self):
        rc = {'int1': self.single_int, 'int2': self.single_int}

        rc_out, = self.dict_of_ints(rc)
        rc_out_2, = self.dict_of_ints(rc_out)
        int1_destructured = rc_out_2['int1']

        with tempfile.TemporaryDirectory() as tempdir:
            in_fp = pathlib.Path(tempdir) / 'int1_destructured.qza'
            int1_destructured.save(in_fp)
            dag = ProvDAG(in_fp)

            out_fp = pathlib.Path(tempdir) / 'rendered.txt'
            out_fn = str(out_fp)
            replay_provenance(ReplayPythonUsage, dag, out_fn)

            with open(out_fp) as fh:
                rendered = fh.read()

        exp = '''\
output_1_artifact_collection, = dummy_plugin_actions.dict_of_ints(
    ints=output_0_artifact_collection,
)
'''
        self.assertIn(exp, rendered)
        self.assertREAppearsOnlyOnce(rendered, r'ResultCollection\(')


class BuildActionUsageTests(CustomAssertions):
    @classmethod
    def setUpClass(cls):
        cls.das = DummyArtifacts()
        cls.tempdir = cls.das.tempdir
        cls.pm = PluginManager()

    @classmethod
    def tearDownClass(cls):
        cls.das.free()

    def test_build_action_usage_python(self):
        plugin = 'dummy_plugin'
        action = 'concatenate_ints'
        cfg = ReplayConfig(
            use=ReplayPythonUsage(),
            use_recorded_metadata=False,
            pm=self.pm
        )

        ns = ReplayNamespaces()
        import_var_1 = ArtifactAPIUsageVariable(
            'imported_ints_0', lambda: None, 'artifact', cfg.use
        )
        import_var_2 = ArtifactAPIUsageVariable(
            'imported_ints_1', lambda: None, 'artifact', cfg.use
        )
        import_uuid_1 = '8dea2f1a-2164-4a85-9f7d-e0641b1db22b'
        import_uuid_2 = '7727c060-5384-445d-b007-b64b41a090ee'
        ns.add_usg_var_record(import_uuid_1, 'imported_ints', import_var_1)
        ns.add_usg_var_record(import_uuid_2, 'imported_ints', import_var_2)

        dag = self.das.concated_ints_v6.dag
        action_uuid = '5035a60e-6f9a-40d4-b412-48ae52255bb5'
        node_uuid = '6facaf61-1676-45eb-ada0-d530be678b27'
        node = dag.get_node_data(node_uuid)
        actions = ActionCollections(
            std_actions={action_uuid: {node_uuid: 'concatenated_ints'}}
        )
        unique_var_name = node.action.output_name + '_0'
        build_action_usage(node, ns, actions.std_actions, action_uuid, cfg)

        usg_var = ns.get_usg_var_record(node_uuid).variable
        out_name = usg_var.to_interface_name()

        self.assertIsInstance(usg_var, UsageVariable)
        self.assertEqual(usg_var.var_type, 'artifact')
        self.assertEqual(usg_var.name, unique_var_name)

        rendered = cfg.use.render()
        self.assertRegex(
            rendered, f"import.*{plugin}.actions as {plugin}_actions"
        )
        self.assertIn(
            f'{out_name}, = dummy_plugin_actions.{action}(', rendered
        )

    def test_build_action_usage_recorded_md(self):
        action = 'identity_with_metadata'
        with tempfile.TemporaryDirectory() as tempdir:
            cfg = ReplayConfig(
                use=ReplayPythonUsage(),
                use_recorded_metadata=False,
                pm=self.pm,
                md_out_dir=tempdir
            )

            action_uuid = '8dae7a81-83ce-48db-9313-6e3131b0933c'
            node_uuid = 'be472b56-d205-43ee-8180-474da575c4d5'

            dag = self.das.concated_ints_with_md.dag
            node = dag.get_node_data(node_uuid)

            ns = ReplayNamespaces()
            mapping_var = ArtifactAPIUsageVariable(
                'imported_mapping_0', lambda: None, 'artifact', cfg.use
            )
            intseq_var_1 = ArtifactAPIUsageVariable(
                'imported_ints_0', lambda: None, 'artifact', cfg.use
            )
            intseq_var_2 = ArtifactAPIUsageVariable(
                'imported_ints_1', lambda: None, 'artifact', cfg.use
            )
            mapping_import_uuid = '8f71b73d-b028-4cbc-9894-738bdfe718bf'
            intseq_import_uuid_1 = '0bb6d731-155a-4dd0-8a1e-98827bc4e0bf'
            intseq_import_uuid_2 = 'e6b37bae-3a14-40f7-87b4-52cf5c7c7a1d'
            ns.add_usg_var_record(
                mapping_import_uuid, 'imported_mapping', mapping_var
            )
            ns.add_usg_var_record(
                intseq_import_uuid_1, 'imported_ints', intseq_var_1
            )
            ns.add_usg_var_record(
                intseq_import_uuid_2, 'imported_ints', intseq_var_2
            )

            actions = ActionCollections(
                std_actions={action_uuid: {node_uuid: 'out'}}
            )
            build_action_usage(node, ns, actions.std_actions, action_uuid, cfg)

            usg_var = ns.get_usg_var_record(node_uuid).variable

            self.assertIsInstance(usg_var, UsageVariable)
            self.assertEqual(usg_var.var_type, 'artifact')
            self.assertEqual(usg_var.name, 'out_0')

            rendered = cfg.use.render()
            self.assertIn('from qiime2 import Metadata', rendered)
            self.assertIn('.view(Metadata)', rendered)
            self.assertIn(f'.{action}(', rendered)


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
    @classmethod
    def setUpClass(cls):
        cls.das = DummyArtifacts()
        cls.tempdir = cls.das.tempdir
        cls.pm = PluginManager()

    @classmethod
    def tearDownClass(cls):
        cls.das.free()

    def test_dedupe_citations(self):
        fn = os.path.join(self.das.datadir, 'dupes.bib')
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
            'ID': 'view|types:2021.2.0|pandas.core.series:Series|0'
        }

        df_20 = {
            'year': ' 2010 ',
            'title': ' Data Structures for Statistical Computing in Python ',
            'pages': ' 51 -- 56 ',
            'editor': ' Stéfan van der Walt and Jarrod Millman ',
            'booktitle': ' Proceedings of the 9th Python in Science Conferen',
            'author': ' Wes McKinney ',
            'ENTRYTYPE': 'inproceedings',
            'ID': 'view|types:2020.2.0|pandas.core.frame:DataFrame|0'
        }

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
                  'method:parse_silva_taxonomy|0'
        }
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
                  'method:parse_silva_taxonomy|1'
        }
        deduped = dedupe_citations([s0, s1])
        self.assertEqual(len(deduped), 2)

    def test_collect_citations_no_dedupe(self):
        dag = self.das.concated_ints_v6.dag
        exp_keys = {
            'framework|qiime2:2023.5.1|0',
            'action|dummy-plugin:0.0.0-dev|method:concatenate_ints|0',
            'plugin|dummy-plugin:0.0.0-dev|0',
            'plugin|dummy-plugin:0.0.0-dev|1',
            'view|dummy-plugin:0.0.0-dev|IntSequenceDirectoryFormat|0',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceDirectoryFormat|0',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceV2DirectoryFormat|0',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceV2DirectoryFormat|1',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceV2DirectoryFormat|2',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceV2DirectoryFormat|3',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceV2DirectoryFormat|4',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceV2DirectoryFormat|5',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceV2DirectoryFormat|6',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceV2DirectoryFormat|7',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceV2DirectoryFormat|8',
        }
        citations = collect_citations(dag, deduplicate=False)
        keys = set(citations.entries_dict.keys())
        self.assertEqual(len(keys), len(exp_keys))
        self.assertEqual(keys, exp_keys)

    def test_collect_citations_dedupe(self):
        dag = self.das.concated_ints_v6.dag
        exp_keys = {
            'framework|qiime2:2023.5.1|0',
            'action|dummy-plugin:0.0.0-dev|method:concatenate_ints|0',
            'plugin|dummy-plugin:0.0.0-dev|0',
            'plugin|dummy-plugin:0.0.0-dev|1',
            'view|dummy-plugin:0.0.0-dev|IntSequenceDirectoryFormat|0',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceDirectoryFormat|0',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceV2DirectoryFormat|4',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceV2DirectoryFormat|5',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceV2DirectoryFormat|6',
            'transformer|dummy-plugin:0.0.0-dev|'
            'builtins:list->IntSequenceV2DirectoryFormat|8'
        }

        citations = collect_citations(dag, deduplicate=True)
        print(citations.entries_dict.keys())
        keys = set(citations.entries_dict.keys())
        self.assertEqual(len(keys), len(exp_keys))
        self.assertEqual(keys, exp_keys)

    def test_collect_citations_no_prov(self):
        dag = self.das.table_v0.dag

        exp_keys = set()
        citations = collect_citations(dag)
        keys = set(citations.entries_dict.keys())
        self.assertEqual(len(keys), 0)
        self.assertEqual(keys, exp_keys)

    def test_replay_citations(self):
        dag = self.das.concated_ints_v6.dag
        exp_keys = {
            'framework|qiime2:2023.5.1|0',
            'action|dummy-plugin:0.0.0-dev|method:concatenate_ints|0',
            'plugin|dummy-plugin:0.0.0-dev|0',
            'plugin|dummy-plugin:0.0.0-dev|1',
            'view|dummy-plugin:0.0.0-dev|IntSequenceDirectoryFormat|0',
        }

        with tempfile.TemporaryDirectory() as tempdir:
            out_fp = os.path.join(tempdir, 'citations.bib')
            replay_citations(dag, out_fp)
            with open(out_fp, 'r') as fp:
                written = fp.read()
                for key in exp_keys:
                    self.assertIn(key, written)

    def test_replay_citations_no_prov(self):
        dag = self.das.table_v0.dag
        exp = "No citations were registered"
        with tempfile.TemporaryDirectory() as tempdir:
            out_fp = os.path.join(tempdir, 'citations.bib')
            replay_citations(dag, out_fp)
            with open(out_fp, 'r') as fp:
                written = fp.read()
                self.assertIn(exp, written)
