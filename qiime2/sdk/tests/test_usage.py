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
        # TODO standardize temporary directories created by QIIME 2
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        self.plugin = get_dummy_plugin()

    def tearDown(self):
        self.test_dir.cleanup()


# TODO
class TestExecutionUsage(TestCaseUsage):
    pass


class TestUsage(TestCaseUsage):
    def test_basic(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.DiagnosticUsage()
        action.examples['concatenate_ints_simple'](use)

        self.assertEqual(2, len(use._recorder))

        obs1, obs2 = use._recorder

        self.assertEqual('comment', obs1['type'])
        self.assertEqual('action', obs2['type'])

        self.assertTrue('basic usage' in obs1['text'])

        self.assertEqual('dummy_plugin', obs2['action'].plugin_id)
        self.assertEqual('concatenate_ints', obs2['action'].action_name)
        self.assertEqual({'int1': 4, 'int2': 2, 'ints1': 'ints_a',
                          'ints2': 'ints_b', 'ints3': 'ints_c'},
                         obs2['input_opts'])
        self.assertEqual({'ints_d': 'concatenated_ints'}, obs2['output_opts'])

    def test_chained(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.DiagnosticUsage()
        action.examples['concatenate_ints_complex'](use)

        self.assertEqual(4, len(use._recorder))

        obs1, obs2, obs3, obs4 = use._recorder

        self.assertEqual('comment', obs1['type'])
        self.assertEqual('action', obs2['type'])
        self.assertEqual('comment', obs3['type'])
        self.assertEqual('action', obs4['type'])

        self.assertTrue('chained usage (pt 1)' in obs1['text'])

        self.assertEqual('dummy_plugin', obs2['action'].plugin_id)
        self.assertEqual('concatenate_ints', obs2['action'].action_name)
        self.assertEqual({'int1': 4, 'int2': 2, 'ints1': 'ints_a',
                          'ints2': 'ints_b', 'ints3': 'ints_c'},
                         obs2['input_opts'])
        self.assertEqual({'ints_d': 'concatenated_ints'}, obs2['output_opts'])

        self.assertTrue('chained usage (pt 2)' in obs3['text'])

        self.assertEqual('dummy_plugin', obs4['action'].plugin_id)
        self.assertEqual('concatenate_ints', obs4['action'].action_name)
        self.assertEqual({'int1': 41, 'int2': 0, 'ints1': 'concatenated_ints',
                          'ints2': 'ints_b', 'ints3': 'ints_c'},
                         obs4['input_opts'])
        self.assertEqual({'concatenated_ints': 'concatenated_ints'},
                         obs4['output_opts'])


class TestUsageAction(TestCaseUsage):
    def test_successful_init(self):
        obs = usage.UsageAction('foo', 'bar')
        self.assertEqual('foo', obs.plugin_id)
        self.assertEqual('bar', obs.action_name)

    def test_invalid_plugin_id(self):
        with self.assertRaisesRegex(ValueError,
                                    'specify a value for plugin_id'):
            usage.UsageAction('', 'bar')

    def test_invalid_action_name(self):
        with self.assertRaisesRegex(ValueError,
                                    'specify a value for action_name'):
            usage.UsageAction('foo', '')

    def test_successful_get_action(self):
        ua = usage.UsageAction('dummy_plugin', 'concatenate_ints')
        obs_action_f, obs_sig = ua.get_action()

        self.assertTrue(isinstance(obs_action_f, action.Method))
        self.assertTrue(isinstance(obs_sig, signature.MethodSignature))

    def test_unknown_action_get_action(self):
        ua = usage.UsageAction('dummy_plugin', 'concatenate_spleens')
        with self.assertRaisesRegex(KeyError,
                                    'No action.*concatenate_spleens'):
            ua.get_action()

    def test_validate_invalid_inputs(self):
        ua = usage.UsageAction('dummy_plugin', 'concatenate_ints')
        with self.assertRaisesRegex(TypeError, 'instance of UsageInputs'):
            ua.validate({}, usage.UsageOutputNames())

    def test_validate_invalid_outputs(self):
        ua = usage.UsageAction('dummy_plugin', 'concatenate_ints')
        with self.assertRaisesRegex(TypeError, 'instance of UsageOutputNames'):
            ua.validate(usage.UsageInputs(), {})


class TestUsageInputs(TestCaseUsage):
    def setUp(self):
        super().setUp()

        def foo(x: dict, z: str, optional: str = None) -> dict:
            return x

        self.signature = signature.MethodSignature(
            foo,
            inputs={'x': Mapping},
            parameters={'z': plugin.Str, 'optional': plugin.Str},
            outputs=[('y', Mapping)],
        )

    def test_successful_init(self):
        obs = usage.UsageInputs(foo='bar')
        self.assertEqual(['foo'], list(obs.values.keys()))
        self.assertEqual(['bar'], list(obs.values.values()))

    def test_validate_missing_input(self):
        ui = usage.UsageInputs(y='hello')

        with self.assertRaisesRegex(ValueError, 'Missing input.*x'):
            ui.validate(self.signature)

    def test_validate_missing_parameter(self):
        ui = usage.UsageInputs(x='hello')

        with self.assertRaisesRegex(ValueError, 'Missing parameter.*z'):
            ui.validate(self.signature)

    def test_validate_extra_values(self):
        ui = usage.UsageInputs(x='hello', z='goodbye', foo=True)

        with self.assertRaisesRegex(ValueError,
                                    'Extra input.*parameter.*foo'):
            ui.validate(self.signature)

    def test_validated_missing_optional_value(self):
        ui = usage.UsageInputs(x='hello', z='goodbye')
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

        with self.assertRaisesRegex(ValueError, 'SDK.*missing output.*hello'):
            uo.validate_derived({'goodbye': 'val'})

    def test_validate_derived_extra_output(self):
        uo = usage.UsageOutputNames(x='goodbye', y='hello')

        with self.assertRaisesRegex(ValueError, 'SDK.*extra output.*peanut'):
            uo.validate_derived({'goodbye': 'val',
                                 'hello': 'val',
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
                                    "Record for 'peanut' not found in scope."):
            use.get_result('peanut')

    def test_action_invalid_action_provided(self):
        use = self.Usage()
        with self.assertRaisesRegex(TypeError, 'provide.*UsageAction'):
            use.action({}, {}, {})
