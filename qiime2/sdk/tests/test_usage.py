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
        records = use.recorder

        self.assertEqual(5, len(records))

        obs1, obs2, obs3, obs4, obs5 = records

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
        records = use.recorder

        self.assertEqual(7, len(records))

        obs1, obs2, obs3, obs4, obs5, obs6, obs7 = records

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
        records = use.recorder

        self.assertEqual(2, len(records))

        obs1, obs2 = records

        self.assertEqual('comment', obs1.source)
        self.assertEqual('comment', obs2.source)

        self.assertEqual('comment 1', obs1.variable)
        self.assertEqual('comment 2', obs2.variable)

    def test_metadata_merging(self):
        action = self.plugin.actions['identity_with_metadata']
        use = usage.DiagnosticUsage()
        action.examples['identity_with_metadata_merging'](use)
        records = use.recorder

        self.assertEqual(5, len(records))

        obs1, obs2, obs3, obs4, obs5 = records

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
        records = use.recorder

        self.assertEqual(4, len(records))

        obs1, obs2, obs3, obs4 = records

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
        records = use.recorder

        self.assertEqual(5, len(records))

        obs1, obs2, obs3, obs4, obs5 = records

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


class TestUsageBaseClass(TestCaseUsage):
    def setUp(self):
        super().setUp()

        class Usage(usage.Usage):
            pass
        self.Usage = Usage

    def test_action_invalid_action_provided(self):
        use = self.Usage()
        with self.assertRaisesRegex(TypeError, 'provide.*UsageAction'):
            use.action({}, {}, {})

    def test_merge_metadata_one_input(self):
        use = self.Usage()
        with self.assertRaisesRegex(ValueError, 'two or more'):
            use.merge_metadata('foo')


class TestUsageVariable(TestCaseUsage):
    def test_basic(self):
        # TODO
        ...


class TestBaseUsage(TestCaseUsage):
    def test_async(self):
        # TODO
        ...


class TestExecutionUsage(TestCaseUsage):
    def test_basic(self):
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

        ints_a, ints_b, single_int1, single_int2, out = use.recorder.values()

        self.assertIsInstance(ints_a.value, Artifact)
        self.assertIsInstance(ints_b.value, Artifact)
        self.assertIsInstance(single_int1.value, Artifact)
        self.assertIsInstance(single_int2.value, Artifact)
        self.assertIsInstance(out.value, Artifact)

    def test_variadic_input_simple_async(self):
        use = usage.ExecutionUsage(asynchronous=True)
        action = self.plugin.actions['variadic_input_method']
        action.examples['variadic_input_simple'](use)

        ints_a, ints_b, single_int1, single_int2, out = use.recorder.values()

        self.assertIsInstance(ints_a.value, Artifact)
        self.assertIsInstance(ints_b.value, Artifact)
        self.assertIsInstance(single_int1.value, Artifact)
        self.assertIsInstance(single_int2.value, Artifact)
        self.assertIsInstance(out.value, Artifact)
