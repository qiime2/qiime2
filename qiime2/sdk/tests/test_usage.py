# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest.mock as mock
import unittest
import tempfile

from qiime2.core.testing.util import get_dummy_plugin
import qiime2.core.testing.examples as examples
from qiime2.sdk import usage, action, UninitializedPluginManagerError
from qiime2 import Metadata, Artifact, MetadataColumn


class TestCaseUsage(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        self.plugin = get_dummy_plugin()

    def tearDown(self):
        self.test_dir.cleanup()


class TestAssertUsageVarType(TestCaseUsage):
    def test_success(self):
        var = usage.UsageVariable('a', lambda: None, 'artifact', None)
        usage.assert_usage_var_type(var, 'artifact')
        self.assertTrue(True)

    def test_failure(self):
        var = usage.UsageVariable('a', lambda: None, 'artifact', None)
        with self.assertRaisesRegex(AssertionError,
                                    'Incorrect.*a,.*visualization.*artifact'):
            usage.assert_usage_var_type(var, 'visualization')


class TestUsageAction(TestCaseUsage):
    def test_successful_init(self):
        obs = usage.UsageAction(plugin_id='foo', action_id='bar')
        self.assertEqual('foo', obs.plugin_id)
        self.assertEqual('bar', obs.action_id)

    def test_invalid_plugin_id(self):
        with self.assertRaisesRegex(ValueError,
                                    'specify a value for plugin_id'):
            usage.UsageAction(plugin_id='', action_id='bar')

    def test_invalid_action_id(self):
        with self.assertRaisesRegex(ValueError,
                                    'specify a value for action_id'):
            usage.UsageAction(plugin_id='foo', action_id='')

    def test_successful_get_action(self):
        ua = usage.UsageAction(
            plugin_id='dummy_plugin', action_id='concatenate_ints')
        obs_action_f = ua.get_action()

        self.assertTrue(isinstance(obs_action_f, action.Method))

    def test_unknown_action_get_action(self):
        ua = usage.UsageAction(
            plugin_id='dummy_plugin', action_id='concatenate_spleens')
        with self.assertRaisesRegex(KeyError,
                                    'No action.*concatenate_spleens'):
            ua.get_action()

    @mock.patch('qiime2.sdk.PluginManager.reuse_existing',
                side_effect=UninitializedPluginManagerError)
    def test_uninitialized_plugin_manager(self, _):
        with self.assertRaisesRegex(UninitializedPluginManagerError,
                                    'create an instance of sdk.PluginManager'):
            usage.UsageAction(
                plugin_id='dummy_plugin', action_id='concatenate_ints')


class TestUsageInputs(TestCaseUsage):
    def test_successful_init(self):
        obs = usage.UsageInputs(foo='bar')
        self.assertEqual(['foo'], list(obs.values.keys()))
        self.assertEqual(['bar'], list(obs.values.values()))


class TestUsageOutputNames(TestCaseUsage):
    def test_successful_init(self):
        obs = usage.UsageOutputNames(foo='bar')
        self.assertEqual(['foo'], list(obs.values.keys()))
        self.assertEqual(['bar'], list(obs.values.values()))

    def test_invalid_init(self):
        with self.assertRaisesRegex(TypeError, 'key.*foo.*string, not.*bool'):
            usage.UsageOutputNames(foo=True)


class TestUsageBaseClass(TestCaseUsage):
    def setUp(self):
        super().setUp()

    def _reset_usage_variables(self, variables):
        for variable in variables:
            variable.value = usage.UsageVariable.DEFERRED

    def test_action_invalid_action_provided(self):
        use = usage.Usage()
        with self.assertRaisesRegex(ValueError, 'expected.*UsageAction'):
            use.action({}, {}, {})

    def test_merge_metadata_one_input(self):
        use = usage.Usage()
        with self.assertRaisesRegex(ValueError, 'two or more'):
            use.merge_metadata('foo')

    def test_action_cache_is_working(self):
        use = usage.Usage()

        ints = use.init_artifact('ints', examples.ints1_factory)
        mapper = use.init_artifact('mapper', examples.mapping1_factory)

        obs = use.action(
            use.UsageAction(plugin_id='dummy_plugin',
                            action_id='typical_pipeline'),
            use.UsageInputs(int_sequence=ints, mapping=mapper,
                            do_extra_thing=True),
            use.UsageOutputNames(out_map='out_map', left='left', right='right',
                                 left_viz='left_viz', right_viz='right_viz')
        )

        # nothing has been executed yet...
        self.assertEqual(obs._cache_info().misses, 0)
        self.assertEqual(obs._cache_info().hits, 0)

        obs_uuids = set()
        for result in obs:
            obs_result = result.execute()
            obs_uuids.add(obs_result.uuid)

        self.assertEqual(len(obs_uuids), 5)

        self.assertEqual(obs._cache_info().misses, 1)
        # 5 results, executed once, minus 1 miss
        self.assertEqual(obs._cache_info().hits, 5 - 1)

        # keep the lru cache intact, but reset the usage variables
        self._reset_usage_variables(obs)

        for result in obs:
            obs_result = result.execute()
            obs_uuids.add(obs_result.uuid)

        # the theory here is that if the memoized action execution wasn't
        # working, we would wind up with twice as many uuids
        self.assertEqual(len(obs_uuids), 5)

        self.assertEqual(obs._cache_info().misses, 1)
        # 5 results, executed twice, minus 1 miss
        self.assertEqual(obs._cache_info().hits, 5 * 2 - 1)

        # this time, reset the lru cache and watch as the results are
        # recompputed
        obs._cache_reset()
        self._reset_usage_variables(obs)

        for result in obs:
            obs_result = result.execute()
            obs_uuids.add(obs_result.uuid)

        # okay, now we should have duplicates of our 5 results
        self.assertEqual(len(obs_uuids), 5 * 2)

        self.assertEqual(obs._cache_info().misses, 1)
        # 5 results, executed once, minus 1 miss
        self.assertEqual(obs._cache_info().hits, 5 - 1)


class TestUsageVariable(TestCaseUsage):
    def test_basic(self):
        # TODO
        ...


class TestDiagnosticUsage(TestCaseUsage):
    def test_basic(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.DiagnosticUsage()
        action.examples['concatenate_ints_simple'](use)

        self.assertEqual(5, len(use.render()))

        obs1, obs2, obs3, obs4, obs5 = use.render()

        self.assertEqual('init_artifact', obs1.source)
        self.assertEqual('init_artifact', obs2.source)
        self.assertEqual('init_artifact', obs3.source)
        self.assertEqual('comment', obs4.source)
        self.assertEqual('action', obs5.source)

        self.assertEqual('ints_a', obs1.variable.name)
        self.assertEqual('ints_b', obs2.variable.name)
        self.assertEqual('ints_c', obs3.variable.name)
        self.assertEqual('This example demonstrates basic usage.',
                         obs4.variable)
        self.assertEqual('ints_d', obs5.variable[0].name)

        self.assertEqual('artifact', obs1.variable.var_type)
        self.assertEqual('artifact', obs2.variable.var_type)
        self.assertEqual('artifact', obs3.variable.var_type)
        self.assertEqual('artifact', obs5.variable[0].var_type)

        self.assertTrue(obs1.variable.is_deferred)
        self.assertTrue(obs2.variable.is_deferred)
        self.assertTrue(obs3.variable.is_deferred)
        self.assertTrue(obs5.variable[0].is_deferred)

    def test_chained(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.DiagnosticUsage()
        action.examples['concatenate_ints_complex'](use)

        self.assertEqual(7, len(use.render()))

        obs1, obs2, obs3, obs4, obs5, obs6, obs7 = use.render()

        self.assertEqual('init_artifact', obs1.source)
        self.assertEqual('init_artifact', obs2.source)
        self.assertEqual('init_artifact', obs3.source)
        self.assertEqual('comment', obs4.source)
        self.assertEqual('action', obs5.source)
        self.assertEqual('comment', obs6.source)
        self.assertEqual('action', obs7.source)

        self.assertEqual('ints_a', obs1.variable.name)
        self.assertEqual('ints_b', obs2.variable.name)
        self.assertEqual('ints_c', obs3.variable.name)
        self.assertEqual('This example demonstrates chained usage (pt 1).',
                         obs4.variable)
        self.assertEqual('ints_d', obs5.variable[0].name)
        self.assertEqual('This example demonstrates chained usage (pt 2).',
                         obs6.variable)
        self.assertEqual('concatenated_ints', obs7.variable[0].name)

        self.assertEqual('artifact', obs1.variable.var_type)
        self.assertEqual('artifact', obs2.variable.var_type)
        self.assertEqual('artifact', obs3.variable.var_type)
        self.assertEqual('artifact', obs5.variable[0].var_type)
        self.assertEqual('artifact', obs7.variable[0].var_type)

        self.assertTrue(obs1.variable.is_deferred)
        self.assertTrue(obs2.variable.is_deferred)
        self.assertTrue(obs3.variable.is_deferred)
        self.assertTrue(obs5.variable[0].is_deferred)
        self.assertTrue(obs7.variable[0].is_deferred)

    def test_comments_only(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.DiagnosticUsage()
        action.examples['comments_only'](use)

        self.assertEqual(2, len(use.render()))

        obs1, obs2 = use.render()

        self.assertEqual('comment', obs1.source)
        self.assertEqual('comment', obs2.source)

        self.assertEqual('comment 1', obs1.variable)
        self.assertEqual('comment 2', obs2.variable)

    def test_metadata_merging(self):
        action = self.plugin.actions['identity_with_metadata']
        use = usage.DiagnosticUsage()
        action.examples['identity_with_metadata_merging'](use)

        self.assertEqual(5, len(use.render()))

        obs1, obs2, obs3, obs4, obs5 = use.render()

        self.assertEqual('init_artifact', obs1.source)
        self.assertEqual('init_metadata', obs2.source)
        self.assertEqual('init_metadata', obs3.source)
        self.assertEqual('merge_metadata', obs4.source)
        self.assertEqual('action', obs5.source)

        self.assertEqual('ints', obs1.variable.name)
        self.assertEqual('md1', obs2.variable.name)
        self.assertEqual('md2', obs3.variable.name)
        self.assertEqual('md3', obs4.variable.name)
        self.assertEqual('out', obs5.variable[0].name)

        self.assertEqual('artifact', obs1.variable.var_type)
        self.assertEqual('metadata', obs2.variable.var_type)
        self.assertEqual('metadata', obs3.variable.var_type)
        self.assertEqual('metadata', obs4.variable.var_type)
        self.assertEqual('artifact', obs5.variable[0].var_type)

        self.assertTrue(obs1.variable.is_deferred)
        self.assertTrue(obs2.variable.is_deferred)
        self.assertTrue(obs3.variable.is_deferred)
        self.assertTrue(obs4.variable.is_deferred)
        self.assertTrue(obs5.variable[0].is_deferred)

    def test_get_metadata_column(self):
        action = self.plugin.actions['identity_with_metadata_column']
        use = usage.DiagnosticUsage()
        action.examples['identity_with_metadata_column_get_mdc'](use)

        self.assertEqual(4, len(use.render()))

        obs1, obs2, obs3, obs4 = use.render()

        self.assertEqual('init_artifact', obs1.source)
        self.assertEqual('init_metadata', obs2.source)
        self.assertEqual('get_metadata_column', obs3.source)
        self.assertEqual('action', obs4.source)

        self.assertEqual('ints', obs1.variable.name)
        self.assertEqual('md', obs2.variable.name)
        self.assertEqual('mdc', obs3.variable.name)
        self.assertEqual('out', obs4.variable[0].name)

        self.assertEqual('artifact', obs1.variable.var_type)
        self.assertEqual('metadata', obs2.variable.var_type)
        self.assertEqual('column', obs3.variable.var_type)
        self.assertEqual('artifact', obs4.variable[0].var_type)

        self.assertTrue(obs1.variable.is_deferred)
        self.assertTrue(obs2.variable.is_deferred)
        self.assertTrue(obs3.variable.is_deferred)
        self.assertTrue(obs4.variable[0].is_deferred)

    def test_optional_inputs(self):
        action = self.plugin.actions['optional_artifacts_method']
        use = usage.DiagnosticUsage()
        action.examples['optional_inputs'](use)

        self.assertEqual(5, len(use.render()))

        obs1, obs2, obs3, obs4, obs5 = use.render()

        self.assertEqual('init_artifact', obs1.source)
        self.assertEqual('action', obs2.source)
        self.assertEqual('action', obs3.source)
        self.assertEqual('action', obs4.source)
        self.assertEqual('action', obs5.source)

        self.assertEqual('ints', obs1.variable.name)
        self.assertEqual('output1', obs2.variable[0].name)
        self.assertEqual('output2', obs3.variable[0].name)
        self.assertEqual('output3', obs4.variable[0].name)
        self.assertEqual('output4', obs5.variable[0].name)

        self.assertEqual('artifact', obs1.variable.var_type)
        self.assertEqual('artifact', obs2.variable[0].var_type)
        self.assertEqual('artifact', obs3.variable[0].var_type)
        self.assertEqual('artifact', obs4.variable[0].var_type)
        self.assertEqual('artifact', obs5.variable[0].var_type)

        self.assertTrue(obs1.variable.is_deferred)
        self.assertTrue(obs2.variable[0].is_deferred)
        self.assertTrue(obs3.variable[0].is_deferred)
        self.assertTrue(obs4.variable[0].is_deferred)
        self.assertTrue(obs5.variable[0].is_deferred)


class TestExecutionUsage(TestCaseUsage):
    def test_basic(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.ExecutionUsage()
        action.examples['concatenate_ints_simple'](use)

        # TODO
        ...

    def test_pipeline(self):
        action = self.plugin.actions['typical_pipeline']
        use = usage.ExecutionUsage()
        action.examples['typical_pipeline_simple'](use)

        # TODO
        ...

    def test_merge_metadata(self):
        use = usage.ExecutionUsage()
        md1 = use.init_metadata('md1', examples.md1_factory)
        md2 = use.init_metadata('md2', examples.md2_factory)
        merged = use.merge_metadata('md3', md1, md2)
        self.assertIsInstance(merged.execute(), Metadata)

    def test_variadic_input_simple(self):
        use = usage.ExecutionUsage()
        action = self.plugin.actions['variadic_input_method']
        action.examples['variadic_input_simple'](use)

        ints_a, ints_b, single_int1, single_int2, out = use.render().values()

        self.assertIsInstance(ints_a.value, Artifact)
        self.assertIsInstance(ints_b.value, Artifact)
        self.assertIsInstance(single_int1.value, Artifact)
        self.assertIsInstance(single_int2.value, Artifact)
        self.assertIsInstance(out.value, Artifact)

    def test_variadic_input_simple_async(self):
        use = usage.ExecutionUsage(asynchronous=True)
        action = self.plugin.actions['variadic_input_method']
        action.examples['variadic_input_simple'](use)

        ints_a, ints_b, single_int1, single_int2, out = use.render().values()

        self.assertIsInstance(ints_a.value, Artifact)
        self.assertIsInstance(ints_b.value, Artifact)
        self.assertIsInstance(single_int1.value, Artifact)
        self.assertIsInstance(single_int2.value, Artifact)
        self.assertIsInstance(out.value, Artifact)

    def test_init_artifact_from_url_error(self):
        use = usage.ExecutionUsage()

        with self.assertRaisesRegex(ValueError, 'Could no.*not-a-url'):
            use.init_artifact_from_url(
                'bad_url_artifact',
                'https://not-a-url.qiime2.org/junk.qza',)

    def test_init_metadata_from_url_error(self):
        use = usage.ExecutionUsage()

        with self.assertRaisesRegex(ValueError, 'Could no.*https://not-a-url'):
            use.init_metadata_from_url(
                'bad_url_metadata',
                'https://not-a-url.qiime2.org/junk.tsv',)

    # def _test_init_artifact_from_url(self):
    #     TODO: need a url to an artifact that the test suite plugin manager
    #     knows about.
    #     artifact_url = ''
    #     use = usage.ExecutionUsage()

    #     a = use.init_artifact_from_url('a', artifact_url)

    #     self.assertIsInstance(a, Artifact)

    def test_init_artifact_from_url_error_on_non_artifact(self):
        # TODO: is this a reliable enough url for tests?
        metadata_url = \
            'https://data.qiime2.org/2022.8/tutorials/' \
            'moving-pictures/sample_metadata.tsv'
        use = usage.ExecutionUsage()

        with self.assertRaisesRegex(ValueError, "Could not.*\n.*a QIIME arc"):
            use.init_artifact_from_url('a', metadata_url)

    def test_init_metadata_from_url_error_on_non_metadata(self):
        url = 'https://www.qiime2.org/'
        use = usage.ExecutionUsage()

        with self.assertRaisesRegex(ValueError, "Could not.*\n.*nized ID"):
            use.init_metadata_from_url('a', url)

    def test_init_metadata_from_url(self):
        # TODO: is this a reliable enough url for tests?
        metadata_url = \
            'https://data.qiime2.org/2022.8/tutorials/' \
            'moving-pictures/sample_metadata.tsv'
        use = usage.ExecutionUsage()

        md = use.init_metadata_from_url('md', metadata_url)

        self.assertIsInstance(md.value, Metadata)

    def test_init_metadata_column_from_url(self):
        # TODO: is this a reliable enough url for tests?
        metadata_url = \
            'https://data.qiime2.org/2022.8/tutorials/' \
            'moving-pictures/sample_metadata.tsv'
        use = usage.ExecutionUsage()

        md = use.init_metadata_from_url('md', metadata_url, column='body-site')

        self.assertIsInstance(md.value, MetadataColumn)

    def test_init_metadata_from_url_epoch(self):
        # TODO: is this a reliable enough url for tests?
        metadata_url = \
            'https://data.qiime2.org/{epoch}/tutorials/' \
            'moving-pictures/sample_metadata.tsv'
        use = usage.ExecutionUsage()

        md = use.init_metadata_from_url('md', metadata_url)

        self.assertIsInstance(md.value, Metadata)

        with self.assertRaisesRegex(ValueError, 'Could no.*{epoch}'):
            use.init_metadata_from_url(
                'bad_url_metadata',
                metadata_url,
                replace_url_epoch=False)
