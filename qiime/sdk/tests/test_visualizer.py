# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import concurrent.futures
import inspect
import os.path
import tempfile
import unittest
import uuid

import qiime.plugin
import qiime.core.type
from qiime.core.callable import Results
from qiime.core.type import VisualizerSignature
from qiime.sdk import Artifact, Visualization, Provenance, Action

from qiime.core.testing.type import IntSequence1, IntSequence2, Mapping
from qiime.core.testing.util import get_dummy_plugin, ArchiveTestingMixin


class TestVisualizer(unittest.TestCase, ArchiveTestingMixin):
    def setUp(self):
        # TODO standardize temporary directories created by QIIME
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        self.plugin = get_dummy_plugin()

    def tearDown(self):
        self.test_dir.cleanup()

    def test_private_constructor(self):
        with self.assertRaisesRegex(NotImplementedError,
                                    'Action constructor.*private'):
            Action()

    def test_from_function_with_artifacts_and_parameters(self):
        visualizer = self.plugin.visualizers['mapping_viz']

        self.assertEqual(visualizer.id, 'mapping_viz')

        exp_sig = VisualizerSignature(
            inputs={
                'mapping1': (Mapping, dict),
                'mapping2': (Mapping, dict)
            },
            parameters={
                'key_label': (qiime.plugin.Str, str),
                'value_label': (qiime.plugin.Str, str)
            },
            defaults={},
            outputs=collections.OrderedDict([
                ('visualization', (qiime.core.type.Visualization, None))
            ])
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
            inputs={
                'ints': (IntSequence1 | IntSequence2, collections.Counter)
            },
            parameters={},
            defaults={},
            outputs=collections.OrderedDict([
                ('visualization', (qiime.core.type.Visualization, None))
            ])
        )
        self.assertEqual(visualizer.signature, exp_sig)

        self.assertEqual(visualizer.name, 'Visualize most common integers')
        self.assertTrue(
            visualizer.description.startswith('This visualizer produces HTML '
                                              'and TSV'))
        self.assertTrue(
            visualizer.source.startswith('\n```python\ndef most_common_viz('))

    def test_is_callable(self):
        self.assertTrue(callable(self.plugin.visualizers['mapping_viz']))
        self.assertTrue(callable(self.plugin.visualizers['most_common_viz']))

    def test_callable_properties(self):
        mapping_viz = self.plugin.visualizers['mapping_viz']
        most_common_viz = self.plugin.visualizers['most_common_viz']

        for visualizer in mapping_viz, most_common_viz:
            self.assertEqual(visualizer.__call__.__name__, '__call__')
            self.assertEqual(visualizer.__call__.__annotations__, {})
            self.assertFalse(hasattr(visualizer.__call__, '__wrapped__'))

    def test_async_properties(self):
        mapping_viz = self.plugin.visualizers['mapping_viz']
        most_common_viz = self.plugin.visualizers['most_common_viz']

        for visualizer in mapping_viz, most_common_viz:
            self.assertEqual(visualizer.async.__name__, 'async')
            self.assertEqual(visualizer.async.__annotations__, {})
            self.assertFalse(hasattr(visualizer.async, '__wrapped__'))

    def test_callable_and_async_signature(self):
        mapping_viz = self.plugin.visualizers['mapping_viz']

        for callable_attr in '__call__', 'async':
            signature = inspect.Signature.from_callable(
                getattr(mapping_viz, callable_attr))
            parameters = list(signature.parameters.items())

            kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
            exp_parameters = [
                ('mapping1', inspect.Parameter('mapping1', kind)),
                ('mapping2', inspect.Parameter('mapping2', kind)),
                ('key_label', inspect.Parameter('key_label', kind)),
                ('value_label', inspect.Parameter('value_label', kind))
            ]
            self.assertEqual(parameters, exp_parameters)

            self.assertEqual(signature.return_annotation,
                             inspect.Signature.empty)

    def test_callable_and_async_different_signature(self):
        # Test that a different Visualizer object has a different dynamic
        # signature.
        most_common_viz = self.plugin.visualizers['most_common_viz']

        for callable_attr in '__call__', 'async':
            signature = inspect.Signature.from_callable(
                getattr(most_common_viz, callable_attr))
            parameters = list(signature.parameters.items())

            kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
            exp_parameters = [
                ('ints', inspect.Parameter('ints', kind))
            ]
            self.assertEqual(parameters, exp_parameters)

            self.assertEqual(signature.return_annotation,
                             inspect.Signature.empty)

    def test_call_with_artifacts_and_parameters(self):
        mapping_viz = self.plugin.visualizers['mapping_viz']

        artifact1 = Artifact._from_view(Mapping, {'foo': 'abc', 'bar': 'def'},
                                        dict, None)
        artifact2 = Artifact._from_view(Mapping, {'baz': 'abc', 'bazz': 'ghi'},
                                        dict, None)

        result = mapping_viz(artifact1, artifact2, 'Key', 'Value')

        # Test properties of the `Results` object.
        self.assertIsInstance(result, tuple)
        self.assertIsInstance(result, Results)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.visualization, result[0])

        result = result[0]

        self.assertIsInstance(result, Visualization)
        self.assertEqual(result.type, qiime.core.type.Visualization)

        provenance = result.provenance
        self.assertIsInstance(provenance, Provenance)
        self.assertIsInstance(provenance.execution_uuid, uuid.UUID)
        self.assertTrue(
            provenance.executor_reference.startswith(mapping_viz.id))
        self.assertEqual(provenance.artifact_uuids, {
            'mapping1': artifact1.uuid,
            'mapping2': artifact2.uuid
        })
        self.assertEqual(provenance.parameter_references, {
            'key_label': 'Key',
            'value_label': 'Value'
        })

        self.assertIsInstance(result.uuid, uuid.UUID)

        # TODO qiime.sdk.Visualization doesn't have an API to access its
        # contents yet. For now, save and assert the correct files are present.
        filepath = os.path.join(self.test_dir.name, 'visualization.qzv')
        result.save(filepath)

        root_dir = str(result.uuid)
        expected = {
            'VERSION',
            'metadata.yaml',
            'README.md',
            'data/index.html',
            'data/css/style.css'
        }

        self.assertArchiveMembers(filepath, root_dir, expected)

    def test_call_with_no_parameters(self):
        most_common_viz = self.plugin.visualizers['most_common_viz']

        artifact = Artifact._from_view(IntSequence1, [42, 42, 10, 0, 42, 5, 0],
                                       list, None)

        result = most_common_viz(artifact)

        # Test properties of the `Results` object.
        self.assertIsInstance(result, tuple)
        self.assertIsInstance(result, Results)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.visualization, result[0])

        result = result[0]

        self.assertIsInstance(result, Visualization)
        self.assertEqual(result.type, qiime.core.type.Visualization)

        provenance = result.provenance
        self.assertIsInstance(provenance, Provenance)
        self.assertIsInstance(provenance.execution_uuid, uuid.UUID)
        self.assertTrue(
            provenance.executor_reference.startswith(most_common_viz.id))
        self.assertEqual(provenance.artifact_uuids, {
            'ints': artifact.uuid
        })
        self.assertEqual(provenance.parameter_references, {})

        self.assertIsInstance(result.uuid, uuid.UUID)

        # TODO qiime.sdk.Visualization doesn't have an API to access its
        # contents yet. For now, save and assert the correct files are present.
        filepath = os.path.join(self.test_dir.name, 'visualization.qzv')
        result.save(filepath)

        root_dir = str(result.uuid)
        expected = {
            'VERSION',
            'metadata.yaml',
            'README.md',
            'data/index.html',
            'data/index.tsv'
        }

        self.assertArchiveMembers(filepath, root_dir, expected)

    def test_async(self):
        mapping_viz = self.plugin.visualizers['mapping_viz']

        artifact1 = Artifact._from_view(Mapping, {'foo': 'abc', 'bar': 'def'},
                                        dict, None)
        artifact2 = Artifact._from_view(Mapping, {'baz': 'abc', 'bazz': 'ghi'},
                                        dict, None)

        future = mapping_viz.async(artifact1, artifact2, 'Key', 'Value')

        self.assertIsInstance(future, concurrent.futures.Future)
        result = future.result()

        # Test properties of the `Results` object.
        self.assertIsInstance(result, tuple)
        self.assertIsInstance(result, Results)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.visualization, result[0])

        result = result[0]

        self.assertIsInstance(result, Visualization)
        self.assertEqual(result.type, qiime.core.type.Visualization)

        provenance = result.provenance
        self.assertIsInstance(provenance, Provenance)
        self.assertIsInstance(provenance.execution_uuid, uuid.UUID)
        self.assertTrue(
            provenance.executor_reference.startswith(mapping_viz.id))
        self.assertEqual(provenance.artifact_uuids, {
            'mapping1': artifact1.uuid,
            'mapping2': artifact2.uuid
        })
        self.assertEqual(provenance.parameter_references, {
            'key_label': 'Key',
            'value_label': 'Value'
        })

        self.assertIsInstance(result.uuid, uuid.UUID)

        # TODO qiime.sdk.Visualization doesn't have an API to access its
        # contents yet. For now, save and assert the correct files are present.
        filepath = os.path.join(self.test_dir.name, 'visualization.qzv')
        result.save(filepath)

        root_dir = str(result.uuid)
        expected = {
            'VERSION',
            'metadata.yaml',
            'README.md',
            'data/index.html',
            'data/css/style.css'
        }

        self.assertArchiveMembers(filepath, root_dir, expected)

    def test_visualizer_callable_output(self):
        artifact = Artifact._from_view(
            Mapping, {'foo': 'abc', 'bar': 'def'}, dict, None)

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


if __name__ == '__main__':
    unittest.main()
