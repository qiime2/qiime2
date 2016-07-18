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
import pkg_resources
import unittest
import uuid

import qiime.plugin
from qiime.core.type import Signature
from qiime.sdk import Artifact, Method, Provenance

from qiime.core.testing.type import IntSequence1, IntSequence2, Mapping
from qiime.core.testing.util import get_dummy_plugin


class TestMethod(unittest.TestCase):
    def setUp(self):
        self.plugin = get_dummy_plugin()

        self.concatenate_ints_sig = Signature(
            inputs={
                'ints1': (IntSequence1 | IntSequence2, list),
                'ints2': (IntSequence1, list),
                'ints3': (IntSequence2, list)
            },
            parameters={
                'int1': (qiime.plugin.Int, int),
                'int2': (qiime.plugin.Int, int)
            },
            outputs=collections.OrderedDict([
                ('concatenated_ints', (IntSequence1, list))
            ])
        )

        self.split_ints_sig = Signature(
            inputs={
                'ints': (IntSequence1, list)
            },
            parameters={},
            outputs=collections.OrderedDict([
                ('left', (IntSequence1, list)),
                ('right', (IntSequence1, list))
            ])
        )

    def test_private_constructor(self):
        with self.assertRaisesRegex(NotImplementedError,
                                    'Method constructor.*private'):
            Method()

    def test_from_markdown_with_artifacts_and_parameters(self):
        method = self.plugin.methods['concatenate_ints_markdown']

        self.assertEqual(method.id, 'concatenate_ints_markdown')
        self.assertEqual(method.signature, self.concatenate_ints_sig)
        self.assertEqual(method.name, 'Concatenate integers')
        self.assertTrue(
            method.description.startswith('This method concatenates integers'))
        self.assertTrue(
            method.source.startswith('## Concatenate some integers'))

    def test_from_markdown_with_multiple_outputs(self):
        method = self.plugin.methods['split_ints_markdown']

        self.assertEqual(method.id, 'split_ints_markdown')
        self.assertEqual(method.signature, self.split_ints_sig)
        self.assertEqual(method.name, 'Split sequence of integers in half')
        self.assertTrue(
            method.description.startswith('This method splits a sequence'))
        self.assertTrue(
            method.source.startswith('### Find midpoint'))

    def test_from_function_with_artifacts_and_parameters(self):
        method = self.plugin.methods['concatenate_ints']

        self.assertEqual(method.id, 'concatenate_ints')
        self.assertEqual(method.signature, self.concatenate_ints_sig)
        self.assertEqual(method.name, 'Concatenate integers')
        self.assertTrue(
            method.description.startswith('This method concatenates integers'))
        self.assertTrue(
            method.source.startswith('\n```python\ndef concatenate_ints('))

    def test_from_function_with_multiple_outputs(self):
        method = self.plugin.methods['split_ints']

        self.assertEqual(method.id, 'split_ints')

        exp_sig = Signature(
            inputs={
                'ints': (IntSequence1, list)
            },
            parameters={},
            outputs=collections.OrderedDict([
                ('left', (IntSequence1, list)),
                ('right', (IntSequence1, list))
            ])
        )
        self.assertEqual(method.signature, exp_sig)

        self.assertEqual(method.name, 'Split sequence of integers in half')
        self.assertTrue(
            method.description.startswith('This method splits a sequence'))
        self.assertTrue(
            method.source.startswith('\n```python\ndef split_ints('))

    def test_from_function_without_parameters(self):
        method = self.plugin.methods['merge_mappings']

        self.assertEqual(method.id, 'merge_mappings')

        exp_sig = Signature(
            inputs={
                'mapping1': (Mapping, dict),
                'mapping2': (Mapping, dict)
            },
            parameters={},
            outputs=collections.OrderedDict([
                ('merged_mapping', (Mapping, dict))
            ])
        )
        self.assertEqual(method.signature, exp_sig)

        self.assertEqual(method.name, 'Merge mappings')
        self.assertTrue(
            method.description.startswith('This method merges two mappings'))
        self.assertTrue(
            method.source.startswith('\n```python\ndef merge_mappings('))

    def test_is_callable(self):
        self.assertTrue(callable(self.plugin.methods['concatenate_ints']))
        self.assertTrue(
            callable(self.plugin.methods['concatenate_ints_markdown']))

    def test_callable_properties(self):
        concatenate_ints = self.plugin.methods['concatenate_ints']
        concatenate_ints_markdown = self.plugin.methods['concatenate_ints']
        merge_mappings = self.plugin.methods['merge_mappings']

        for method in (concatenate_ints, concatenate_ints_markdown,
                       merge_mappings):
            self.assertEqual(method.__call__.__name__, '__call__')
            self.assertEqual(method.__call__.__annotations__, {})
            self.assertFalse(hasattr(method.__call__, '__wrapped__'))

    def test_async_properties(self):
        concatenate_ints = self.plugin.methods['concatenate_ints']
        concatenate_ints_markdown = self.plugin.methods['concatenate_ints']
        merge_mappings = self.plugin.methods['merge_mappings']

        for method in (concatenate_ints, concatenate_ints_markdown,
                       merge_mappings):
            self.assertEqual(method.async.__name__, 'async')
            self.assertEqual(method.async.__annotations__, {})
            self.assertFalse(hasattr(method.async, '__wrapped__'))

    def test_callable_and_async_signature_with_artifacts_and_parameters(self):
        # Signature with input artifacts and parameters (i.e. primitives).
        concatenate_ints = self.plugin.methods['concatenate_ints']
        concatenate_ints_markdown = self.plugin.methods['concatenate_ints']

        for method in concatenate_ints, concatenate_ints_markdown:
            for callable_attr in '__call__', 'async':
                signature = inspect.Signature.from_callable(
                    getattr(method, callable_attr))
                parameters = list(signature.parameters.items())

                kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
                exp_parameters = [
                    ('ints1', inspect.Parameter('ints1', kind)),
                    ('ints2', inspect.Parameter('ints2', kind)),
                    ('ints3', inspect.Parameter('ints3', kind)),
                    ('int1', inspect.Parameter('int1', kind)),
                    ('int2', inspect.Parameter('int2', kind))
                ]
                self.assertEqual(parameters, exp_parameters)

                self.assertEqual(signature.return_annotation,
                                 inspect.Signature.empty)

    def test_callable_and_async_signature_with_no_parameters(self):
        # Signature without parameters (i.e. primitives), only input artifacts.
        method = self.plugin.methods['merge_mappings']

        for callable_attr in '__call__', 'async':
            signature = inspect.Signature.from_callable(
                getattr(method, callable_attr))
            parameters = list(signature.parameters.items())

            kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
            exp_parameters = [
                ('mapping1', inspect.Parameter('mapping1', kind)),
                ('mapping2', inspect.Parameter('mapping2', kind))
            ]
            self.assertEqual(parameters, exp_parameters)

            self.assertEqual(signature.return_annotation,
                             inspect.Signature.empty)

    def test_call_with_artifacts_and_parameters(self):
        concatenate_ints = self.plugin.methods['concatenate_ints']
        concatenate_ints_markdown = \
            self.plugin.methods['concatenate_ints_markdown']

        artifact1 = Artifact._from_view([0, 42, 43], IntSequence1, None)
        artifact2 = Artifact._from_view([99, -22], IntSequence2, None)

        for method in concatenate_ints, concatenate_ints_markdown:
            result = method(artifact1, artifact1, artifact2, 55, 1)

            self.assertIsInstance(result, Artifact)
            self.assertEqual(result.type, IntSequence1)

            provenance = result.provenance
            self.assertIsInstance(provenance, Provenance)
            self.assertIsInstance(provenance.execution_uuid, uuid.UUID)
            self.assertTrue(
                provenance.executor_reference.startswith(method.id))
            self.assertEqual(provenance.artifact_uuids, {
                'ints1': artifact1.uuid,
                'ints2': artifact1.uuid,
                'ints3': artifact2.uuid
            })
            self.assertEqual(provenance.parameter_references, {
                'int1': '55',
                'int2': '1'
            })

            self.assertIsInstance(result.uuid, uuid.UUID)

            # Can retrieve multiple views of different type.
            exp_list_view = [0, 42, 43, 0, 42, 43, 99, -22, 55, 1]
            self.assertEqual(result.view(list), exp_list_view)
            self.assertEqual(result.view(list), exp_list_view)

            exp_counter_view = collections.Counter(
                {0: 2, 42: 2, 43: 2, 99: 1, -22: 1, 55: 1, 1: 1})
            self.assertEqual(result.view(collections.Counter),
                             exp_counter_view)
            self.assertEqual(result.view(collections.Counter),
                             exp_counter_view)

            # Accepts IntSequence1 | IntSequence2
            artifact3 = Artifact._from_view([10, 20], IntSequence2, None)
            result = method(artifact3, artifact1, artifact2, 55, 1)

            self.assertEqual(result.type, IntSequence1)
            self.assertEqual(result.view(list),
                             [10, 20, 0, 42, 43, 99, -22, 55, 1])

    def test_call_with_multiple_outputs(self):
        split_ints = self.plugin.methods['split_ints']
        split_ints_markdown = self.plugin.methods['split_ints_markdown']

        artifact = Artifact._from_view([0, 42, -2, 43, 6], IntSequence1, None)

        for method in split_ints, split_ints_markdown:
            result = method(artifact)

            self.assertIsInstance(result, tuple)
            self.assertEqual(len(result), 2)

            for output_artifact in result:
                self.assertEqual(output_artifact.type, IntSequence1)

                provenance = output_artifact.provenance
                self.assertIsInstance(provenance, Provenance)
                self.assertIsInstance(provenance.execution_uuid, uuid.UUID)
                self.assertTrue(
                    provenance.executor_reference.startswith(method.id))
                self.assertEqual(provenance.artifact_uuids, {
                    'ints': artifact.uuid
                })
                self.assertEqual(provenance.parameter_references, {})

                self.assertIsInstance(output_artifact.uuid, uuid.UUID)

            # Output artifacts have the same provenance.
            self.assertEqual(result[0].provenance, result[1].provenance)

            # Output artifacts have different UUIDs.
            self.assertNotEqual(result[0].uuid, result[1].uuid)

            self.assertEqual(result[0].view(list), [0, 42])
            self.assertEqual(result[1].view(list), [-2, 43, 6])

    def test_call_with_no_parameters(self):
        merge_mappings = self.plugin.methods['merge_mappings']

        artifact1 = Artifact._from_view({'foo': 'abc', 'bar': 'def'}, Mapping,
                                        None)
        artifact2 = Artifact._from_view({'bazz': 'abc'}, Mapping, None)

        result = merge_mappings(artifact1, artifact2)

        self.assertIsInstance(result, Artifact)
        self.assertEqual(result.type, Mapping)

        provenance = result.provenance
        self.assertIsInstance(provenance, Provenance)
        self.assertIsInstance(provenance.execution_uuid, uuid.UUID)
        self.assertTrue(
            provenance.executor_reference.startswith(merge_mappings.id))
        self.assertEqual(provenance.artifact_uuids, {
            'mapping1': artifact1.uuid,
            'mapping2': artifact2.uuid
        })
        self.assertEqual(provenance.parameter_references, {})

        self.assertIsInstance(result.uuid, uuid.UUID)

        self.assertEqual(result.view(dict),
                         {'foo': 'abc', 'bar': 'def', 'bazz': 'abc'})

    def test_async(self):
        concatenate_ints = self.plugin.methods['concatenate_ints']
        concatenate_ints_markdown = \
            self.plugin.methods['concatenate_ints_markdown']

        artifact1 = Artifact._from_view([0, 42, 43], IntSequence1, None)
        artifact2 = Artifact._from_view([99, -22], IntSequence2, None)

        for method in concatenate_ints, concatenate_ints_markdown:
            future = method.async(artifact1, artifact1, artifact2, 55, 1)

            self.assertIsInstance(future, concurrent.futures.Future)
            result = future.result()

            self.assertIsInstance(result, Artifact)
            self.assertEqual(result.type, IntSequence1)

            provenance = result.provenance
            self.assertIsInstance(provenance, Provenance)
            self.assertIsInstance(provenance.execution_uuid, uuid.UUID)
            self.assertTrue(
                provenance.executor_reference.startswith(method.id))
            self.assertEqual(provenance.artifact_uuids, {
                'ints1': artifact1.uuid,
                'ints2': artifact1.uuid,
                'ints3': artifact2.uuid
            })
            self.assertEqual(provenance.parameter_references, {
                'int1': '55',
                'int2': '1'
            })

            self.assertIsInstance(result.uuid, uuid.UUID)

            # Can retrieve multiple views of different type.
            exp_list_view = [0, 42, 43, 0, 42, 43, 99, -22, 55, 1]
            self.assertEqual(result.view(list), exp_list_view)
            self.assertEqual(result.view(list), exp_list_view)

            exp_counter_view = collections.Counter(
                {0: 2, 42: 2, 43: 2, 99: 1, -22: 1, 55: 1, 1: 1})
            self.assertEqual(result.view(collections.Counter),
                             exp_counter_view)
            self.assertEqual(result.view(collections.Counter),
                             exp_counter_view)

            # Accepts IntSequence1 | IntSequence2
            artifact3 = Artifact._from_view([10, 20], IntSequence2, None)
            future = method.async(artifact3, artifact1, artifact2, 55, 1)
            result = future.result()

            self.assertEqual(result.type, IntSequence1)
            self.assertEqual(result.view(list),
                             [10, 20, 0, 42, 43, 99, -22, 55, 1])

    def test_markdown_input_section_validation(self):
        fp = pkg_resources.resource_filename(
            'qiime.sdk.tests', 'data/method_bad_input_section.md')

        with self.assertRaisesRegex(TypeError, 'input section'):
            Method.from_markdown(fp)

    def test_markdown_parameters_section_validation(self):
        fp = pkg_resources.resource_filename(
            'qiime.sdk.tests', 'data/method_bad_parameters_section.md')

        with self.assertRaisesRegex(TypeError, 'parameters section'):
            Method.from_markdown(fp)

    def test_markdown_outputs_section_validation(self):
        fp = pkg_resources.resource_filename(
            'qiime.sdk.tests', 'data/method_bad_outputs_section.md')

        with self.assertRaisesRegex(TypeError, 'outputs section'):
            Method.from_markdown(fp)


if __name__ == '__main__':
    unittest.main()
