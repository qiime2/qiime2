# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
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
from qiime2.sdk import usage, action
from qiime2 import plugin


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

        self.assertEqual(5, len(use.recorder))

        obs1, obs2, obs3, obs4, obs5 = use.recorder

        self.assertEqual('init_data', obs1['type'])
        self.assertEqual('init_data', obs2['type'])
        self.assertEqual('init_data', obs3['type'])
        self.assertEqual('comment', obs4['type'])
        self.assertEqual('action', obs5['type'])

        self.assertTrue('basic usage' in obs4['text'])

        self.assertEqual('dummy_plugin', obs5['action'].plugin_id)
        self.assertEqual('concatenate_ints', obs5['action'].action_id)
        self.assertEqual({'int1': 4, 'int2': 2, 'ints1': 'ints_a',
                          'ints2': 'ints_b', 'ints3': 'ints_c'},
                         obs5['input_opts'])
        self.assertEqual({'concatenated_ints': 'ints_d'}, obs5['output_opts'])

    def test_chained(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.DiagnosticUsage()
        action.examples['concatenate_ints_complex'](use)

        self.assertEqual(7, len(use.recorder))

        obs1, obs2, obs3, obs4, obs5, obs6, obs7 = use.recorder

        self.assertEqual('init_data', obs1['type'])
        self.assertEqual('init_data', obs2['type'])
        self.assertEqual('init_data', obs3['type'])
        self.assertEqual('comment', obs4['type'])
        self.assertEqual('action', obs5['type'])
        self.assertEqual('comment', obs6['type'])
        self.assertEqual('action', obs7['type'])

        self.assertTrue('chained usage (pt 1)' in obs4['text'])

        self.assertEqual('dummy_plugin', obs5['action'].plugin_id)
        self.assertEqual('concatenate_ints', obs5['action'].action_id)
        self.assertEqual({'int1': 4, 'int2': 2, 'ints1': 'ints_a',
                          'ints2': 'ints_b', 'ints3': 'ints_c'},
                         obs5['input_opts'])
        self.assertEqual({'concatenated_ints': 'ints_d'}, obs5['output_opts'])

        self.assertTrue('chained usage (pt 2)' in obs6['text'])

        self.assertEqual('dummy_plugin', obs7['action'].plugin_id)
        self.assertEqual('concatenate_ints', obs7['action'].action_id)
        self.assertEqual({'int1': 41, 'int2': 0, 'ints1': 'ints_d',
                          'ints2': 'ints_b', 'ints3': 'ints_c'},
                         obs7['input_opts'])
        self.assertEqual({'concatenated_ints': 'concatenated_ints'},
                         obs7['output_opts'])

    def test_comments_only(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.DiagnosticUsage()
        action.examples['comments_only'](use)

        self.assertEqual(2, len(use.recorder))

        obs1, obs2 = use.recorder

        self.assertEqual('comment', obs1['type'])
        self.assertEqual('comment', obs2['type'])

        self.assertEqual('comment 1', obs1['text'])
        self.assertEqual('comment 2', obs2['text'])

    def test_metadata_merging(self):
        action = self.plugin.actions['identity_with_metadata']
        use = usage.DiagnosticUsage()
        action.examples['identity_with_metadata_merging'](use)

        self.assertEqual(5, len(use.recorder))

        obs1, obs2, obs3, obs4, obs5 = use.recorder

        self.assertEqual('init_data', obs1['type'])
        self.assertEqual('init_data', obs2['type'])
        self.assertEqual('init_data', obs3['type'])
        self.assertEqual('merge_metadata', obs4['type'])
        self.assertEqual('action', obs5['type'])

    def test_get_metadata_column(self):
        action = self.plugin.actions['identity_with_metadata_column']
        use = usage.DiagnosticUsage()
        action.examples['identity_with_metadata_column_get_mdc'](use)

        self.assertEqual(4, len(use.recorder))

        obs1, obs2, obs3, obs4 = use.recorder

        self.assertEqual('init_data', obs1['type'])
        self.assertEqual('init_data', obs2['type'])
        self.assertEqual('get_metadata_column', obs3['type'])
        self.assertEqual('action', obs4['type'])

    def test_use_metadata_column(self):
        action = self.plugin.actions['identity_with_metadata_column']
        use = usage.DiagnosticUsage()
        action.examples['identity_with_metadata_column_from_factory'](use)

        self.assertEqual(3, len(use.recorder))

        obs1, obs2, obs3 = use.recorder

        self.assertEqual('init_data', obs1['type'])
        self.assertEqual('init_data', obs2['type'])
        self.assertEqual('action', obs3['type'])

    def test_optional_inputs(self):
        action = self.plugin.actions['optional_artifacts_method']
        use = usage.DiagnosticUsage()

        action.examples['optional_inputs'](use)

        self.assertEqual(5, len(use.recorder))

        obs1, obs2, obs3, obs4, obs5 = use.recorder
        self.assertEqual('init_data', obs1['type'])
        self.assertEqual('action', obs2['type'])
        self.assertEqual('action', obs3['type'])
        self.assertEqual('action', obs4['type'])
        self.assertEqual('action', obs5['type'])


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
            usage.ScopeRecord('foo', assert_has_line_matching='spleen')
