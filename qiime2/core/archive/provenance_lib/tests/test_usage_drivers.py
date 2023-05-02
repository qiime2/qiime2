import os
import pathlib
import tempfile
import unittest

from qiime2.sdk.usage import UsageAction

from ..replay import replay_provenance
from .test_parse import DATA_DIR
from .._usage_drivers import (
    MissingPluginError, ReplayCLIUsage, _get_action_if_plugin_present,
)


class MiscHelperFunctionTests(unittest.TestCase):
    def test_get_action_if_plugin_present_plugin_present(self):
        real_action = UsageAction('diversity', 'core_metrics')
        action = _get_action_if_plugin_present(real_action)
        self.assertEqual('diversity', action.plugin_id)
        self.assertEqual('core_metrics', action.id)

    def test_get_action_if_plugin_present_plugin_missing(self):
        fake_action = UsageAction('imaginary', 'action')
        with self.assertRaisesRegex(
                MissingPluginError,
                "(?s)missing one or more plugins.*library"):
            _get_action_if_plugin_present(fake_action)


class ReplayCLIUsageTests(unittest.TestCase):
    def test_init_metadata(self):
        use = ReplayCLIUsage()
        var = use.init_metadata(name='testing', factory=lambda: None)
        print(var)
        self.assertEqual(var.name, '<your metadata filepath>')
        self.assertEqual(var.var_type, 'metadata')

    def test_init_metadata_with_dumped_md_fn(self):
        use = ReplayCLIUsage()
        var = use.init_metadata(
            name='testing', factory=lambda: None, dumped_md_fn='some_md')
        self.assertEqual(var.var_type, 'metadata')
        self.assertEqual(var.name, '"some_md.tsv"')


class ReplayPythonUsageTests(unittest.TestCase):
    def test_template_action_lumps_many_outputs(self):
        """
        ReplayPythonUsage._template_action should "lump" multiple outputs from
        one command into a single Results-like object when the total number of
        outputs from a single command > 5

        In these cases, our rendering should look like:
        `action_results = plugin_actions.action()...`
        instead of:
        `_, _, thing3, _, _, _ = plugin_actions.action()...`

        In this artifact, we are only replaying one results from core-metrics,
        but because core_metrics has a million results it should stil lump em.
        """
        in_fp = os.path.join(DATA_DIR, 'v5_uu_emperor.qzv')
        driver = 'python3'
        exp = ('(?s)action_results = diversity_actions.core_metrics_phylo.*'
               'unweighted_unifrac_emperor.*action_results.unweighted_unifrac')
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = pathlib.Path(tmpdir) / 'action_collection.txt'
            replay_provenance(in_fp, out_path, driver)

            with open(out_path, 'r') as fp:
                rendered = fp.read()
        self.assertRegex(rendered, exp)

    def test_template_action_does_not_lump_four_outputs(self):
        """
        ReplayPythonUsage._template_action should not "lump" multiple outputs
        one command into a single Results-like object when the total number of
        outputs from a single command <= 5, unless the total number of results
        is high (see above).

        In these cases, our rendering should look like:
        `_, _, thing3, _ = plugin_actions.action()...`
        instead of:
        `action_results = plugin_actions.action()...`

        In this case, we are replaying one result from an action which has four
        results. It should not lump em.
        """
        in_fp = os.path.join(DATA_DIR, 'v5_uu_emperor.qzv')
        driver = 'python3'
        exp = ('(?s)_, _, _, rooted_tree_0 = phylogeny_actions.align_to_tre.*')
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = pathlib.Path(tmpdir) / 'action_collection.txt'
            replay_provenance(in_fp, out_path, driver)

            with open(out_path, 'r') as fp:
                rendered = fp.read()
        self.assertRegex(rendered, exp)

    def test_template_action_lumps_three_variables(self):
        """
        ReplayPythonUsage._template_action should "lump" multiple outputs from
        one command into a single Results-like object when there are more than
        two usage variables (i.e. replay of 3+ results from a single command)

        In these cases, our rendering should look like:
        ```
        action_results = plugin_actions.action(...)
        thing1 = action_results.thinga
        etc.
        ```
        instead of:
        `thing1, _, thing3, _, thing5, _ = plugin_actions.action()...`

        In this test, we are replaying three results from dada2.denoise_single,
        which should be lumped.
        """
        in_fp = os.path.join(DATA_DIR, 'lump_three_vars_test')
        driver = 'python3'
        e1 = ('action_results = dada2_actions.denoise_single')
        e2 = ('representative_sequences_0 = action_results.representative_seq')
        e3 = ('denoising_stats_0 = action_results.denoising_stats')
        e4 = ('table_0 = action_results.table')
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = pathlib.Path(tmpdir) / 'action_collection.txt'
            replay_provenance(in_fp, out_path, driver)

            with open(out_path, 'r') as fp:
                rendered = fp.read()
        self.assertRegex(rendered, e1)
        self.assertRegex(rendered, e2)
        self.assertRegex(rendered, e3)
        self.assertRegex(rendered, e4)

    def test_template_action_does_not_lump_two_vars(self):
        """
        ReplayPythonUsage._template_action should not "lump" multiple outputs
        from one command into a single Results-like object when the total count
        of usage variables (i.e. replayed outputs) from a single command < 3,
        unless the total number of outputs is high (see above).

        In these cases, our rendering should look like:
        `thing1, _, thing3, _ = plugin_actions.action()...`
        instead of:
        `action_results = plugin_actions.action()...`

        In this case, we are replaying two results from dada2.denoise_single,
        which should not be lumped.
        """
        in_fp = os.path.join(DATA_DIR, 'v5_uu_emperor.qzv')
        driver = 'python3'
        exp = ('(?s)table_0, representative_sequences_0, _ = dada2_actions.*')
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = pathlib.Path(tmpdir) / 'action_collection.txt'
            replay_provenance(in_fp, out_path, driver)

            with open(out_path, 'r') as fp:
                rendered = fp.read()
        self.assertRegex(rendered, exp)


class ActionPatchTests(unittest.TestCase):
    def test_missing_plugin(self):
        """
        action_patch raises a MissingPluginError if a plugin from provenance
        is missing in the Env. The test .qza requires rescript, which is not
        included in the test env.
        """
        in_fp = os.path.join(DATA_DIR, 'rescript-based-taxonomy.qza')
        exp = ('(?s)QIIME 2 deployment.*missing.*plugins.*rescript')
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = pathlib.Path(tmpdir) / 'whatever.thing'
            with self.assertRaisesRegex(MissingPluginError, exp):
                replay_provenance(in_fp, out_path)
