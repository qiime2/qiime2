# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import inspect

import pandas as pd

import qiime2
import qiime2.sdk
from qiime2.core.testing.util import get_dummy_plugin
from qiime2.core.testing.type import IntSequence1, SingleInt, Mapping
from qiime2.plugin import Visualization, Int, Bool


class TestPipeline(unittest.TestCase):
    def setUp(self):
        self.plugin = get_dummy_plugin()
        self.single_int = qiime2.Artifact.import_data(SingleInt, -1)
        self.int_sequence = qiime2.Artifact.import_data(IntSequence1,
                                                        [1, 2, 3])
        self.mapping = qiime2.Artifact.import_data(Mapping, {'foo': '42'})

    def test_private_constructor(self):
        with self.assertRaisesRegex(NotImplementedError,
                                    'Pipeline constructor.*private'):
            qiime2.sdk.Pipeline()

    def test_from_function_spot_check(self):
        typical_pipeline = self.plugin.pipelines['typical_pipeline']
        self.assertEqual(typical_pipeline.id, 'typical_pipeline')

        assert typical_pipeline.signature.inputs
        for spec in typical_pipeline.signature.inputs.values():
            assert spec.has_description()
            assert spec.has_qiime_type()
            assert not spec.has_view_type()
            assert not spec.has_default()

        spec = typical_pipeline.signature.parameters['add']
        assert spec.has_default()

    def test_from_function_optional(self):
        optional_artifact_pipeline = self.plugin.pipelines[
            'optional_artifact_pipeline']

        spec = optional_artifact_pipeline.signature.inputs['single_int']
        assert spec.has_default()

    def test_is_callable(self):
        assert callable(self.plugin.pipelines['typical_pipeline'])

    def test_callable_and_async_signature(self):
        # Shouldn't include `ctx`
        typical_pipeline = self.plugin.pipelines['typical_pipeline']
        kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
        exp_parameters = [
            ('int_sequence', inspect.Parameter(
                'int_sequence', kind, annotation=IntSequence1)),
            ('mapping', inspect.Parameter(
                'mapping', kind, annotation=Mapping)),
            ('do_extra_thing', inspect.Parameter(
                'do_extra_thing', kind, annotation=Bool)),
            ('add', inspect.Parameter(
                'add', kind, default=1, annotation=Int))
        ]

        for callable_attr in '__call__', 'asynchronous':
            signature = inspect.Signature.from_callable(
                getattr(typical_pipeline, callable_attr))
            parameters = list(signature.parameters.items())

            self.assertEqual(parameters, exp_parameters)

    def test_signatures_independent(self):
        typical_pipeline = self.plugin.pipelines['typical_pipeline']
        parameter_only_pipeline = self.plugin.pipelines[
            'parameter_only_pipeline']

        for callable_attr in '__call__', 'asynchronous':
            signature_a = inspect.Signature.from_callable(
                getattr(typical_pipeline, callable_attr))

            signature_b = inspect.Signature.from_callable(
                getattr(parameter_only_pipeline, callable_attr))

            self.assertNotEqual(signature_a, signature_b)

    def iter_callables(self, name):
        pipeline = self.plugin.pipelines[name]
        yield pipeline
        yield lambda *args, **kwargs: pipeline.asynchronous(
            *args, **kwargs).result()

    def test_parameter_only_pipeline(self):
        index = pd.Index(['a', 'b', 'c'], name='id', dtype=object)
        df = pd.DataFrame({'col1': ['2', '1', '3']}, index=index, dtype=object)
        metadata = qiime2.Metadata(df)
        for call in self.iter_callables('parameter_only_pipeline'):
            results = call(100)
            self.assertEqual(results.foo.view(list), [100, 2, 3])
            self.assertEqual(results.bar.view(list),
                             [100, 2, 3, 100, 2, 3, 100, 2, 3, 100, 2])

            results = call(3, int2=4, metadata=metadata)
            self.assertEqual(results.foo.view(list), [3, 4, 3])
            self.assertEqual(results.bar.view(list),
                             [3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4])

    def test_typical_pipeline(self):
        for call in self.iter_callables('typical_pipeline'):
            results = call(self.int_sequence, self.mapping, False)

            self.assertEqual(results.left_viz.type, Visualization)
            self.assertEqual(results.left.view(list), [1])
            self.assertEqual(results.right.view(list), [2, 3])
            self.assertNotEqual(results.out_map.uuid, self.mapping.uuid)
            self.assertEqual(results.out_map.view(dict),
                             self.mapping.view(dict))

            results = call(self.int_sequence, self.mapping, True, add=5)
            self.assertEqual(results.left.view(list), [6])
            self.assertEqual(results.right.view(list), [2, 3])

            with self.assertRaisesRegex(ValueError, 'Bad mapping'):
                m = qiime2.Artifact.import_data(Mapping, {'a': 1})
                call(self.int_sequence, m, False)

    def test_optional_artifact_pipeline(self):
        for call in self.iter_callables('optional_artifact_pipeline'):
            ints, = call(self.int_sequence)
            self.assertEqual(ints.view(list), [1, 2, 3, 4])

            ints, = call(self.int_sequence, single_int=self.single_int)
            self.assertEqual(ints.view(list), [1, 2, 3, -1])

    def test_visualizer_only_pipeline(self):
        for call in self.iter_callables('visualizer_only_pipeline'):
            viz1, viz2 = call(self.mapping)

            self.assertEqual(viz1.type, Visualization)
            self.assertEqual(viz2.type, Visualization)

    def test_pipeline_in_pipeline(self):
        for call in self.iter_callables('pipelines_in_pipeline'):
            results = call(self.int_sequence, self.mapping)

            self.assertEqual(results.int1.view(int), 4)
            self.assertEqual(results.right_viz.type, Visualization)
            self.assertEqual(len(results), 8)

            with self.assertRaisesRegex(ValueError, 'Bad mapping'):
                m = qiime2.Artifact.import_data(Mapping, {1: 1})
                call(self.int_sequence, m)

    def test_pointless_pipeline(self):
        for call in self.iter_callables('pointless_pipeline'):
            single_int, = call()
            self.assertEqual(single_int.type, SingleInt)
            self.assertEqual(single_int.view(int), 4)

    def test_failing_from_arity(self):
        for call in self.iter_callables('failing_pipeline'):
            with self.assertRaisesRegex(TypeError, 'match number.*3.*1'):
                call(self.int_sequence, break_from='arity')

    def test_failing_from_return_view(self):
        for call in self.iter_callables('failing_pipeline'):
            with self.assertRaisesRegex(TypeError, 'Result objects'):
                call(self.int_sequence, break_from='return-view')

    def test_failing_from_method(self):
        for call in self.iter_callables('failing_pipeline'):
            with self.assertRaisesRegex(ValueError, "Key 'foo' exists"):
                call(self.int_sequence, break_from='method')

    def test_failing_from_type(self):
        for call in self.iter_callables('failing_pipeline'):
            with self.assertRaisesRegex(TypeError, 'Mapping.*SingleInt'):
                call(self.int_sequence, break_from='type')

    def test_failing_from_internal(self):
        for call in self.iter_callables('failing_pipeline'):
            with self.assertRaisesRegex(ValueError, 'this never works'):
                call(self.int_sequence, break_from='internal')

    def test_failing_from_missing_plugin(self):
        for call in self.iter_callables('failing_pipeline'):
            with self.assertRaisesRegex(ValueError, r'plugin.*not\%a\$plugin'):
                call(self.int_sequence, break_from='no-plugin')

    def test_failing_from_missing_action(self):
        for call in self.iter_callables('failing_pipeline'):
            with self.assertRaisesRegex(ValueError, r'action.*not\%a\$method'):
                call(self.int_sequence, break_from='no-action')


if __name__ == '__main__':
    unittest.main()
