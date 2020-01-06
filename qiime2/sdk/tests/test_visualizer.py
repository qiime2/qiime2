# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import concurrent.futures
import inspect
import os.path
import tempfile
import unittest
import uuid

import qiime2.plugin
import qiime2.core.type
from qiime2.core.type import VisualizerSignature, Str, Range
from qiime2.core.type.visualization import Visualization as VisualizationType
from qiime2.sdk import Artifact, Visualization, Visualizer, Results

from qiime2.core.testing.visualizer import (most_common_viz, mapping_viz,
                                            params_only_viz, no_input_viz)
from qiime2.core.testing.type import IntSequence1, IntSequence2, Mapping
from qiime2.core.testing.util import get_dummy_plugin, ArchiveTestingMixin


class TestVisualizer(unittest.TestCase, ArchiveTestingMixin):
    def setUp(self):
        # TODO standardize temporary directories created by QIIME 2
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        self.plugin = get_dummy_plugin()

    def tearDown(self):
        self.test_dir.cleanup()

    def test_private_constructor(self):
        with self.assertRaisesRegex(NotImplementedError,
                                    'Visualizer constructor.*private'):
            Visualizer()

    def test_from_function_with_artifacts_and_parameters(self):
        visualizer = self.plugin.visualizers['mapping_viz']

        self.assertEqual(visualizer.id, 'mapping_viz')

        exp_sig = VisualizerSignature(
            mapping_viz,
            inputs={
                'mapping1': Mapping,
                'mapping2': Mapping
            },
            parameters={
                'key_label': qiime2.plugin.Str,
                'value_label': qiime2.plugin.Str
            },
        )
        self.assertEqual(visualizer.signature, exp_sig)

        self.assertEqual(visualizer.name, 'Visualize two mappings')
        self.assertTrue(
            visualizer.description.startswith('This visualizer produces an '
                                              'HTML visualization'))
        self.assertTrue(
            visualizer.source.startswith('\n```python\ndef mapping_viz('))

    def test_from_function_without_parameters(self):
        visualizer = self.plugin.visualizers['most_common_viz']

        self.assertEqual(visualizer.id, 'most_common_viz')

        exp_sig = VisualizerSignature(
            most_common_viz,
            inputs={
                'ints': IntSequence1 | IntSequence2
            },
            parameters={}
        )
        self.assertEqual(visualizer.signature, exp_sig)

        self.assertEqual(visualizer.name, 'Visualize most common integers')
        self.assertTrue(
            visualizer.description.startswith('This visualizer produces HTML '
                                              'and TSV'))
        self.assertTrue(
            visualizer.source.startswith('\n```python\ndef most_common_viz('))

    def test_from_function_with_parameters_only(self):
        visualizer = self.plugin.visualizers['params_only_viz']

        self.assertEqual(visualizer.id, 'params_only_viz')

        exp_sig = VisualizerSignature(
            params_only_viz,
            inputs={},
            parameters={
                'name': qiime2.plugin.Str,
                'age': qiime2.plugin.Int % Range(0, None)
            }
        )
        self.assertEqual(visualizer.signature, exp_sig)

        self.assertEqual(visualizer.name, 'Parameters only viz')
        self.assertTrue(
            visualizer.description.startswith('This visualizer only accepts '
                                              'parameters.'))
        self.assertTrue(
            visualizer.source.startswith('\n```python\ndef params_only_viz('))

    def test_from_function_without_inputs_or_parameters(self):
        visualizer = self.plugin.visualizers['no_input_viz']

        self.assertEqual(visualizer.id, 'no_input_viz')

        exp_sig = VisualizerSignature(
            no_input_viz,
            inputs={},
            parameters={}
        )
        self.assertEqual(visualizer.signature, exp_sig)

        self.assertEqual(visualizer.name, 'No input viz')
        self.assertTrue(
            visualizer.description.startswith('This visualizer does not '
                                              'accept any'))
        self.assertTrue(
            visualizer.source.startswith('\n```python\ndef no_input_viz('))

    def test_is_callable(self):
        self.assertTrue(callable(self.plugin.visualizers['mapping_viz']))
        self.assertTrue(callable(self.plugin.visualizers['most_common_viz']))

    def test_callable_properties(self):
        mapping_viz = self.plugin.visualizers['mapping_viz']
        most_common_viz = self.plugin.visualizers['most_common_viz']

        mapping_exp = {
            'mapping1': Mapping, 'return': (VisualizationType,),
            'key_label': Str, 'mapping2': Mapping, 'value_label': Str}
        most_common_exp = {
            'ints': IntSequence1 | IntSequence2,
            'return': (VisualizationType,)}

        mapper = {
            mapping_viz: mapping_exp,
            most_common_viz: most_common_exp}

        for visualizer, exp in mapper.items():
            self.assertEqual(visualizer.__call__.__name__, '__call__')
            self.assertEqual(visualizer.__call__.__annotations__, exp)
            self.assertFalse(hasattr(visualizer.__call__, '__wrapped__'))

    def test_async_properties(self):
        mapping_viz = self.plugin.visualizers['mapping_viz']
        most_common_viz = self.plugin.visualizers['most_common_viz']

        mapping_exp = {
            'mapping1': Mapping, 'return': (VisualizationType,),
            'key_label': Str, 'mapping2': Mapping, 'value_label': Str}
        most_common_exp = {
            'ints': IntSequence1 | IntSequence2,
            'return': (VisualizationType,)}

        mapper = {
            mapping_viz: mapping_exp,
            most_common_viz: most_common_exp}

        for visualizer, exp in mapper.items():
            self.assertEqual(visualizer.asynchronous.__name__, 'asynchronous')
            self.assertEqual(visualizer.asynchronous.__annotations__, exp)
            self.assertFalse(hasattr(visualizer.asynchronous, '__wrapped__'))

    def test_callable_and_async_signature(self):
        mapping_viz = self.plugin.visualizers['mapping_viz']

        for callable_attr in '__call__', 'asynchronous':
            signature = inspect.Signature.from_callable(
                getattr(mapping_viz, callable_attr))
            parameters = list(signature.parameters.items())

            kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
            exp_parameters = [
                ('mapping1', inspect.Parameter(
                    'mapping1', kind, annotation=Mapping)),
                ('mapping2', inspect.Parameter(
                    'mapping2', kind, annotation=Mapping)),
                ('key_label', inspect.Parameter(
                    'key_label', kind, annotation=Str)),
                ('value_label', inspect.Parameter(
                    'value_label', kind, annotation=Str))
            ]

            self.assertEqual(parameters, exp_parameters)

    def test_callable_and_async_different_signature(self):
        # Test that a different Visualizer object has a different dynamic
        # signature.
        most_common_viz = self.plugin.visualizers['most_common_viz']

        for callable_attr in '__call__', 'asynchronous':
            signature = inspect.Signature.from_callable(
                getattr(most_common_viz, callable_attr))
            parameters = list(signature.parameters.items())

            kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
            exp_parameters = [
                ('ints', inspect.Parameter(
                    'ints', kind, annotation=IntSequence1 | IntSequence2))
            ]

            self.assertEqual(parameters, exp_parameters)

    def test_call_with_artifacts_and_parameters(self):
        mapping_viz = self.plugin.visualizers['mapping_viz']

        artifact1 = Artifact.import_data(Mapping, {'foo': 'abc', 'bar': 'def'})
        artifact2 = Artifact.import_data(
            Mapping, {'baz': 'abc', 'bazz': 'ghi'})

        result = mapping_viz(artifact1, artifact2, 'Key', 'Value')

        # Test properties of the `Results` object.
        self.assertIsInstance(result, tuple)
        self.assertIsInstance(result, Results)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.visualization, result[0])

        result = result[0]

        self.assertIsInstance(result, Visualization)
        self.assertEqual(result.type, qiime2.core.type.Visualization)

        self.assertIsInstance(result.uuid, uuid.UUID)

        # TODO qiime2.sdk.Visualization doesn't have an API to access its
        # contents yet. For now, save and assert the correct files are present.
        filepath = os.path.join(self.test_dir.name, 'visualization.qzv')
        result.save(filepath)

        root_dir = str(result.uuid)
        expected = {
            'VERSION',
            'checksums.md5',
            'metadata.yaml',
            'data/index.html',
            'data/css/style.css',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/citations.bib',
            'provenance/action/action.yaml',
            'provenance/artifacts/%s/metadata.yaml' % artifact1.uuid,
            'provenance/artifacts/%s/VERSION' % artifact1.uuid,
            'provenance/artifacts/%s/citations.bib' % artifact1.uuid,
            'provenance/artifacts/%s/action/action.yaml' % artifact1.uuid,
            'provenance/artifacts/%s/metadata.yaml' % artifact2.uuid,
            'provenance/artifacts/%s/VERSION' % artifact2.uuid,
            'provenance/artifacts/%s/citations.bib' % artifact2.uuid,
            'provenance/artifacts/%s/action/action.yaml' % artifact2.uuid
        }

        self.assertArchiveMembers(filepath, root_dir, expected)

    def test_call_with_no_parameters(self):
        most_common_viz = self.plugin.visualizers['most_common_viz']

        artifact = Artifact.import_data(
            IntSequence1, [42, 42, 10, 0, 42, 5, 0])

        result = most_common_viz(artifact)

        # Test properties of the `Results` object.
        self.assertIsInstance(result, tuple)
        self.assertIsInstance(result, Results)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.visualization, result[0])

        result = result[0]

        self.assertIsInstance(result, Visualization)
        self.assertEqual(result.type, qiime2.core.type.Visualization)

        self.assertIsInstance(result.uuid, uuid.UUID)

        # TODO qiime2.sdk.Visualization doesn't have an API to access its
        # contents yet. For now, save and assert the correct files are present.
        filepath = os.path.join(self.test_dir.name, 'visualization.qzv')
        result.save(filepath)

        root_dir = str(result.uuid)
        expected = {
            'VERSION',
            'checksums.md5',
            'metadata.yaml',
            'data/index.html',
            'data/index.tsv',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/citations.bib',
            'provenance/action/action.yaml',
            'provenance/artifacts/%s/metadata.yaml' % artifact.uuid,
            'provenance/artifacts/%s/VERSION' % artifact.uuid,
            'provenance/artifacts/%s/citations.bib' % artifact.uuid,
            'provenance/artifacts/%s/action/action.yaml' % artifact.uuid
        }

        self.assertArchiveMembers(filepath, root_dir, expected)

    def test_call_with_parameters_only(self):
        params_only_viz = self.plugin.visualizers['params_only_viz']

        # Parameters all have default values.
        result, = params_only_viz()

        self.assertIsInstance(result, Visualization)
        self.assertEqual(result.type, qiime2.core.type.Visualization)
        self.assertIsInstance(result.uuid, uuid.UUID)

        filepath = os.path.join(self.test_dir.name, 'visualization.qzv')
        result.save(filepath)

        root_dir = str(result.uuid)
        expected = {
            'VERSION',
            'checksums.md5',
            'metadata.yaml',
            'data/index.html',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/citations.bib',
            'provenance/action/action.yaml'
        }

        self.assertArchiveMembers(filepath, root_dir, expected)

    def test_call_without_inputs_or_parameters(self):
        no_input_viz = self.plugin.visualizers['no_input_viz']

        result, = no_input_viz()

        self.assertIsInstance(result, Visualization)
        self.assertEqual(result.type, qiime2.core.type.Visualization)
        self.assertIsInstance(result.uuid, uuid.UUID)

        filepath = os.path.join(self.test_dir.name, 'visualization.qzv')
        result.save(filepath)

        root_dir = str(result.uuid)
        expected = {
            'VERSION',
            'checksums.md5',
            'metadata.yaml',
            'data/index.html',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/citations.bib',
            'provenance/action/action.yaml'
        }

        self.assertArchiveMembers(filepath, root_dir, expected)

    def test_asynchronous(self):
        mapping_viz = self.plugin.visualizers['mapping_viz']

        artifact1 = Artifact.import_data(Mapping, {'foo': 'abc', 'bar': 'def'})
        artifact2 = Artifact.import_data(
            Mapping, {'baz': 'abc', 'bazz': 'ghi'})

        future = mapping_viz.asynchronous(artifact1, artifact2, 'Key', 'Value')

        self.assertIsInstance(future, concurrent.futures.Future)
        result = future.result()

        # Test properties of the `Results` object.
        self.assertIsInstance(result, tuple)
        self.assertIsInstance(result, Results)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.visualization, result[0])

        result = result[0]

        self.assertIsInstance(result, Visualization)
        self.assertEqual(result.type, qiime2.core.type.Visualization)

        self.assertIsInstance(result.uuid, uuid.UUID)

        # TODO qiime2.sdk.Visualization doesn't have an API to access its
        # contents yet. For now, save and assert the correct files are present.
        filepath = os.path.join(self.test_dir.name, 'visualization.qzv')
        result.save(filepath)

        root_dir = str(result.uuid)
        expected = {
            'VERSION',
            'checksums.md5',
            'metadata.yaml',
            'data/index.html',
            'data/css/style.css',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/citations.bib',
            'provenance/action/action.yaml',
            'provenance/artifacts/%s/metadata.yaml' % artifact1.uuid,
            'provenance/artifacts/%s/VERSION' % artifact1.uuid,
            'provenance/artifacts/%s/citations.bib' % artifact1.uuid,
            'provenance/artifacts/%s/action/action.yaml' % artifact1.uuid,
            'provenance/artifacts/%s/metadata.yaml' % artifact2.uuid,
            'provenance/artifacts/%s/VERSION' % artifact2.uuid,
            'provenance/artifacts/%s/citations.bib' % artifact2.uuid,
            'provenance/artifacts/%s/action/action.yaml' % artifact2.uuid
        }

        self.assertArchiveMembers(filepath, root_dir, expected)

    def test_visualizer_callable_output(self):
        artifact = Artifact.import_data(Mapping, {'foo': 'abc', 'bar': 'def'})

        # Callable returns a value from `return_vals`
        return_vals = (True, False, [], {}, '', 0, 0.0)
        for return_val in return_vals:
            def func(output_dir: str, foo: dict) -> None:
                return return_val

            self.plugin.visualizers.register_function(
                func, {'foo': Mapping}, {}, '', ''
            )
            visualizer = self.plugin.visualizers['func']

            with self.assertRaisesRegex(TypeError, "should not return"):
                visualizer(foo=artifact)

        # Callable returns None (default function return)
        def func(output_dir: str, foo: dict) -> None:
            return None

        self.plugin.visualizers.register_function(
            func, {'foo': Mapping}, {}, '', ''
        )
        visualizer = self.plugin.visualizers['func']

        # Should not raise an exception
        output = visualizer(foo=artifact)
        self.assertIsInstance(output, Results)
        self.assertIsInstance(output.visualization, Visualization)

    def test_docstring(self):
        mapping_viz = self.plugin.visualizers['mapping_viz']
        common_viz = self.plugin.visualizers['most_common_viz']
        params_only_viz = self.plugin.visualizers['params_only_viz']
        no_input_viz = self.plugin.visualizers['no_input_viz']

        obs = mapping_viz.__call__.__doc__
        self.assertEqual(obs, exp_mapping_viz)

        obs = common_viz.__call__.__doc__
        self.assertEqual(obs, exp_common_viz)

        obs = params_only_viz.__call__.__doc__
        self.assertEqual(obs, exp_params_only_viz)

        obs = no_input_viz.__call__.__doc__
        self.assertEqual(obs, exp_no_input_viz)


exp_mapping_viz = """\
Visualize two mappings

This visualizer produces an HTML visualization of two key-value mappings,
each sorted in alphabetical order by key.

Parameters
----------
mapping1 : Mapping
mapping2 : Mapping
key_label : Str
value_label : Str

Returns
-------
visualization : Visualization
"""

exp_common_viz = """\
Visualize most common integers

This visualizer produces HTML and TSV outputs containing the input sequence
of integers ordered from most- to least-frequently occurring, along with
their respective frequencies.

Parameters
----------
ints : IntSequence1 | IntSequence2

Returns
-------
visualization : Visualization
"""

exp_params_only_viz = """\
Parameters only viz

This visualizer only accepts parameters.

Parameters
----------
name : Str, optional
age : Int % Range(0, None), optional

Returns
-------
visualization : Visualization
"""

exp_no_input_viz = """\
No input viz

This visualizer does not accept any type of input.

Returns
-------
visualization : Visualization
"""
if __name__ == '__main__':
    unittest.main()
