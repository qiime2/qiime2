# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import tempfile

from qiime2.core.type import signature
from qiime2.core.testing.util import get_dummy_plugin
from qiime2.core.testing.type import Mapping
import qiime2.core.testing.examples as examples
from qiime2.sdk import usage, action
from qiime2 import plugin, Metadata, Artifact


class TestCaseUsage(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        self.plugin = get_dummy_plugin()

    def tearDown(self):
        self.test_dir.cleanup()


class TestUsage(TestCaseUsage):
    def test_basic(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.DiagnosticUsage()
        action.examples['concatenate_ints_simple'](use)
        records = use._get_records()

        self.assertEqual(5, len(records))

        obs1, obs2, obs3, obs4, obs5 = records.values()

        self.assertEqual('init_data', obs1.source)
        self.assertEqual('init_data', obs2.source)
        self.assertEqual('init_data', obs3.source)
        self.assertEqual('comment', obs4.source)
        self.assertEqual('action', obs5.source)

        self.assertTrue('basic usage' in obs4.result['text'])

        self.assertEqual('dummy_plugin', obs5.result['plugin_id'])
        self.assertEqual('concatenate_ints', obs5.result['action_id'])
        self.assertEqual({'int1': 4, 'int2': 2,
                          'ints1': {'ref': 'ints_a', 'source': 'init_data'},
                          'ints2': {'ref': 'ints_b', 'source': 'init_data'},
                          'ints3': {'ref': 'ints_c', 'source': 'init_data'}},
                         obs5.result['input_opts'])
        self.assertEqual({'concatenated_ints': 'ints_d'},
                         obs5.result['output_opts'])

    def test_chained(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.DiagnosticUsage()
        action.examples['concatenate_ints_complex'](use)
        records = use._get_records()

        self.assertEqual(7, len(records))

        obs1, obs2, obs3, obs4, obs5, obs6, obs7 = records.values()

        self.assertEqual('init_data', obs1.source)
        self.assertEqual('init_data', obs2.source)
        self.assertEqual('init_data', obs3.source)
        self.assertEqual('comment', obs4.source)
        self.assertEqual('action', obs5.source)
        self.assertEqual('comment', obs6.source)
        self.assertEqual('action', obs7.source)

        self.assertTrue('chained usage (pt 1)' in obs4.result['text'])

        self.assertEqual('dummy_plugin', obs5.result['plugin_id'])
        self.assertEqual('concatenate_ints', obs5.result['action_id'])
        self.assertEqual({'int1': 4, 'int2': 2,
                          'ints1': {'ref': 'ints_a', 'source': 'init_data'},
                          'ints2': {'ref': 'ints_b', 'source': 'init_data'},
                          'ints3': {'ref': 'ints_c', 'source': 'init_data'}},
                         obs5.result['input_opts'])
        self.assertEqual({'concatenated_ints': 'ints_d'},
                         obs5.result['output_opts'])

        self.assertTrue('chained usage (pt 2)' in obs6.result['text'])

        self.assertEqual('dummy_plugin', obs7.result['plugin_id'])
        self.assertEqual('concatenate_ints', obs7.result['action_id'])
        exp7 = {'int1': 41, 'int2': 0,
                'ints1': {'action_id': 'concatenate_ints',
                          'input_opts': {'int1': 4, 'int2': 2,
                                         'ints1': {'ref': 'ints_a',
                                                   'source': 'init_data'},
                                         'ints2': {'ref': 'ints_b',
                                                   'source': 'init_data'},
                                         'ints3': {'ref': 'ints_c',
                                                   'source': 'init_data'}},
                          'output_opt': 'concatenated_ints',
                          'output_opts': {'concatenated_ints': 'ints_d'},
                          'plugin_id': 'dummy_plugin',
                          'source': 'action'},
                'ints2': {'ref': 'ints_b', 'source': 'init_data'},
                'ints3': {'ref': 'ints_c', 'source': 'init_data'}}
        self.assertEqual(exp7, obs7.result['input_opts'])
        self.assertEqual({'concatenated_ints': 'concatenated_ints'},
                         obs7.result['output_opts'])

    def test_comments_only(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.DiagnosticUsage()
        action.examples['comments_only'](use)
        records = use._get_records()

        self.assertEqual(2, len(records))

        obs1, obs2 = records.values()

        self.assertEqual('comment', obs1.source)
        self.assertEqual('comment', obs2.source)

        self.assertEqual('comment 1', obs1.result['text'])
        self.assertEqual('comment 2', obs2.result['text'])

    def test_metadata_merging(self):
        action = self.plugin.actions['identity_with_metadata']
        use = usage.DiagnosticUsage()
        action.examples['identity_with_metadata_merging'](use)
        records = use._get_records()

        self.assertEqual(5, len(records))

        obs1, obs2, obs3, obs4, obs5 = records.values()

        self.assertEqual('init_data', obs1.source)
        self.assertEqual('init_metadata', obs2.source)
        self.assertEqual('init_metadata', obs3.source)
        self.assertEqual('merge_metadata', obs4.source)
        self.assertEqual('action', obs5.source)

    def test_get_metadata_column(self):
        action = self.plugin.actions['identity_with_metadata_column']
        use = usage.DiagnosticUsage()
        action.examples['identity_with_metadata_column_get_mdc'](use)
        records = use._get_records()

        self.assertEqual(4, len(records))

        obs1, obs2, obs3, obs4 = records.values()

        self.assertEqual('init_data', obs1.source)
        self.assertEqual('init_metadata', obs2.source)
        self.assertEqual('get_metadata_column', obs3.source)
        self.assertEqual('action', obs4.source)

    def test_use_init_collection_data(self):
        action = self.plugin.actions['variadic_input_method']
        use = usage.DiagnosticUsage()
        action.examples['variadic_input_simple'](use)
        records = use._get_records()

        self.assertEqual(7, len(records))

        obs1, obs2, obs3, obs4, obs5, obs6, obs7 = records.values()

        self.assertEqual('init_data', obs1.source)
        self.assertEqual('init_data', obs2.source)
        self.assertEqual('init_data_collection', obs3.source)

        self.assertEqual('init_data', obs4.source)
        self.assertEqual('init_data', obs5.source)
        self.assertEqual('init_data_collection', obs6.source)
        self.assertEqual('action', obs7.source)

        self.assertEqual(set, type(obs7.result['input_opts']['nums']))

        self.assertEqual('ints', obs7.result['input_opts']['ints'][0]['ref'])
        self.assertEqual('int_set',
                         obs7.result['input_opts']['int_set'][0]['ref'])

    def test_optional_inputs(self):
        action = self.plugin.actions['optional_artifacts_method']
        use = usage.DiagnosticUsage()
        records = use._get_records()

        action.examples['optional_inputs'](use)

        self.assertEqual(3, len(records))

        obs1, obs2, obs3 = records.values()

        self.assertEqual('init_data', obs1.source)
        self.assertEqual('action', obs2.source)
        self.assertEqual('action', obs3.source)


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
        obs_action_f, obs_sig = ua.get_action()

        self.assertTrue(isinstance(obs_action_f, action.Method))
        self.assertTrue(isinstance(obs_sig, signature.MethodSignature))

    def test_unknown_action_get_action(self):
        ua = usage.UsageAction(
            plugin_id='dummy_plugin', action_id='concatenate_spleens')
        with self.assertRaisesRegex(KeyError,
                                    'No action.*concatenate_spleens'):
            ua.get_action()

    def test_validate_invalid_inputs(self):
        ua = usage.UsageAction(
            plugin_id='dummy_plugin', action_id='concatenate_ints')
        with self.assertRaisesRegex(TypeError, 'instance of UsageInputs'):
            ua.validate({}, usage.UsageOutputNames())

    def test_validate_invalid_outputs(self):
        ua = usage.UsageAction(
            plugin_id='dummy_plugin', action_id='concatenate_ints')
        with self.assertRaisesRegex(TypeError, 'instance of UsageOutputNames'):
            ua.validate(usage.UsageInputs(), {})


class TestUsageInputs(TestCaseUsage):
    def setUp(self):
        super().setUp()

        def foo(x: dict, z: str,
                optional_input: dict = None,
                optional_param: str = None) -> dict:
            return x

        self.signature = signature.MethodSignature(
            foo,
            inputs={'x': Mapping, 'optional_input': Mapping},
            parameters={'z': plugin.Str, 'optional_param': plugin.Str},
            outputs=[('y', Mapping)],
        )

    def test_successful_init(self):
        obs = usage.UsageInputs(foo='bar')
        self.assertEqual(['foo'], list(obs.values.keys()))
        self.assertEqual(['bar'], list(obs.values.values()))

    def test_validate_missing_required_input(self):
        ui = usage.UsageInputs(y='hello',
                               optional_input='a', optional_param='b')

        with self.assertRaisesRegex(ValueError, 'Missing input.*x'):
            ui.validate(self.signature)

    def test_validate_missing_required_parameter(self):
        ui = usage.UsageInputs(x='hello',
                               optional_input='a', optional_param='b')

        with self.assertRaisesRegex(ValueError, 'Missing parameter.*z'):
            ui.validate(self.signature)

    def test_validate_extra_values(self):
        ui = usage.UsageInputs(x='hello', z='goodbye', foo=True,
                               optional_input='a', optional_param='b')

        with self.assertRaisesRegex(ValueError,
                                    'Extra input.*parameter.*foo'):
            ui.validate(self.signature)

    def test_validate_missing_optional_input(self):
        ui = usage.UsageInputs(x='hello', z='goodbye', optional_param='a')
        ui.validate(self.signature)
        self.assertTrue(True)

    def test_validate_missing_optional_parameter(self):
        ui = usage.UsageInputs(x='hello', z='goodbye', optional_input='a')
        ui.validate(self.signature)
        self.assertTrue(True)

    def test_type_of_input(self):
        test_inputs = usage.UsageInputs(x=[1, 2, 3], z={7, 8, 9})
        test_inputs.validate(self.signature)

        self.assertIsInstance(test_inputs.values['x'], list)
        self.assertIsInstance(test_inputs.values['z'], set)


class TestUsageOutputNames(TestCaseUsage):
    def setUp(self):
        super().setUp()

        def foo(x: dict, z: str, optional: str = None) -> (dict, dict):
            return x

        self.signature = signature.MethodSignature(
            foo,
            inputs={'x': Mapping},
            parameters={'z': plugin.Str, 'optional': plugin.Str},
            outputs=[('y', Mapping), ('a', Mapping)],
        )

    def test_successful_init(self):
        obs = usage.UsageOutputNames(foo='bar')
        self.assertEqual(['foo'], list(obs.values.keys()))
        self.assertEqual(['bar'], list(obs.values.values()))

    def test_invalid_init(self):
        with self.assertRaisesRegex(TypeError, 'key.*foo.*string, not.*bool'):
            usage.UsageOutputNames(foo=True)

    def test_validate_missing_output(self):
        uo = usage.UsageOutputNames(y='hello')

        with self.assertRaisesRegex(ValueError, 'Missing output.*a'):
            uo.validate(self.signature)

    def test_validate_extra_output(self):
        uo = usage.UsageOutputNames(y='goodbye', a='hello', peanut='noeyes')

        with self.assertRaisesRegex(ValueError, 'Extra output.*peanut'):
            uo.validate(self.signature)

    def test_validate_derived_missing_output(self):
        uo = usage.UsageOutputNames(x='goodbye', y='hello')

        with self.assertRaisesRegex(ValueError, 'SDK.*missing output.*y'):
            uo.validate_computed({'x': 'val'})

    def test_validate_derived_extra_output(self):
        uo = usage.UsageOutputNames(x='goodbye', y='hello')

        with self.assertRaisesRegex(ValueError, 'SDK.*extra output.*peanut'):
            uo.validate_computed({'x': 'val',
                                  'y': 'val',
                                  'peanut': 'val'})


class TestUsageBaseClass(TestCaseUsage):
    def setUp(self):
        super().setUp()

        class Usage(usage.Usage):
            pass
        self.Usage = Usage

    def test_get_result_invalid(self):
        use = self.Usage()
        with self.assertRaisesRegex(KeyError,
                                    'No record with ref id: "peanut"'):
            use.get_result('peanut')

    def test_action_invalid_action_provided(self):
        use = self.Usage()
        with self.assertRaisesRegex(TypeError, 'provide.*UsageAction'):
            use.action({}, {}, {})

    def test_merge_metadata_one_input(self):
        use = self.Usage()
        with self.assertRaisesRegex(ValueError, 'two or more'):
            use.merge_metadata('foo')


class TestScopeRecord(TestCaseUsage):
    def test_invalid_assert_has_line_matching(self):
        with self.assertRaisesRegex(TypeError, 'should be a `callable`'):
            usage.ScopeRecord('foo', 'value', 'source',
                              assert_has_line_matching='spleen')


class TestExecutionUsage(TestCaseUsage):
    def test_init_data(self):
        use = usage.ExecutionUsage()

        with self.assertRaisesRegex(ValueError, 'expected an Artifact'):
            use.init_data('name', lambda: object)

        with self.assertRaisesRegex(ValueError, 'not all .* Artifacts'):
            use.init_data('name', lambda: [object])

        with self.assertRaisesRegex(TypeError, 'expected Metadata'):
            use.init_metadata('name', lambda: object)

        with self.assertRaisesRegex(ValueError, 'expected a ScopeRecord.'):
            use.init_data_collection('', list, object)

        with self.assertRaisesRegex(ValueError, 'expected a ScopeRecord.'):
            use.init_data_collection('', list,
                                     usage.ScopeRecord('', object, ''), object)

    def test_merge_metadata(self):
        use = usage.ExecutionUsage()
        md1 = use.init_metadata('md1', examples.md1_factory)
        md2 = use.init_metadata('md2', examples.md2_factory)
        merged = use.merge_metadata('md3', md1, md2)
        self.assertIsInstance(merged.result, Metadata)

    def test_variadic_input_simple(self):
        use = usage.ExecutionUsage()
        action = self.plugin.actions['variadic_input_method']
        action.examples['variadic_input_simple'](use)
        ints_a = use._get_record('ints_a')
        ints_b = use._get_record('ints_b')
        ints = use._get_record('ints')
        single_int1 = use._get_record('single_int1')
        single_int2 = use._get_record('single_int2')
        int_set = use._get_record('int_set')
        out = use._get_record('out')
        self.assertIsInstance(ints_a.result, Artifact)
        self.assertIsInstance(ints_b.result, Artifact)
        self.assertIsInstance(ints.result, list)
        self.assertEqual(ints.result[0], ints_a.result)
        self.assertEqual(ints.result[1], ints_b.result)
        self.assertIsInstance(single_int1.result, Artifact)
        self.assertIsInstance(single_int2.result, Artifact)
        self.assertIsInstance(int_set.result, set)
        self.assertIn(single_int1.result, int_set.result)
        self.assertIn(single_int2.result, int_set.result)
        self.assertIsInstance(out.result, Artifact)

    def test_variadic_input_simple_async(self):
        use = usage.ExecutionUsage(asynchronous=True)
        action = self.plugin.actions['variadic_input_method']
        action.examples['variadic_input_simple'](use)
        ints_a = use._get_record('ints_a')
        ints_b = use._get_record('ints_b')
        ints = use._get_record('ints')
        single_int1 = use._get_record('single_int1')
        single_int2 = use._get_record('single_int2')
        int_set = use._get_record('int_set')
        out = use._get_record('out')
        self.assertIsInstance(ints_a.result, Artifact)
        self.assertIsInstance(ints_b.result, Artifact)
        self.assertIsInstance(ints.result, list)
        self.assertEqual(ints.result[0], ints_a.result)
        self.assertEqual(ints.result[1], ints_b.result)
        self.assertIsInstance(single_int1.result, Artifact)
        self.assertIsInstance(single_int2.result, Artifact)
        self.assertIsInstance(int_set.result, set)
        self.assertIn(single_int1.result, int_set.result)
        self.assertIn(single_int2.result, int_set.result)
        self.assertIsInstance(out.result, Artifact)
