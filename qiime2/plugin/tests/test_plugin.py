# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import types
import unittest

import qiime2.plugin
import qiime2.sdk

from qiime2.core.testing.type import (IntSequence1, IntSequence2, Mapping,
                                      FourInts, Kennel, Dog, Cat, SingleInt)
from qiime2.core.testing.util import get_dummy_plugin


class TestPlugin(unittest.TestCase):
    def setUp(self):
        self.plugin = get_dummy_plugin()

    def test_name(self):
        self.assertEqual(self.plugin.name, 'dummy-plugin')

    def test_version(self):
        self.assertEqual(self.plugin.version, '0.0.0-dev')

    def test_website(self):
        self.assertEqual(self.plugin.website,
                         'https://github.com/qiime2/qiime2')

    def test_package(self):
        self.assertEqual(self.plugin.package, 'qiime2.core.testing')

    def test_citations(self):
        self.assertEqual(self.plugin.citations[0].type, 'article')

    def test_user_support_text(self):
        self.assertEqual(self.plugin.user_support_text,
                         'For help, see https://qiime2.org')

    def test_short_description_text(self):
        self.assertEqual(self.plugin.short_description,
                         'Dummy plugin for testing.')

    def test_description_text(self):
        self.assertEqual(self.plugin.description,
                         'Description of dummy plugin.')

    def test_citations_default(self):
        plugin = qiime2.plugin.Plugin(
            name='local-dummy-plugin',
            version='0.0.0-dev',
            website='https://github.com/qiime2/qiime2',
            package='qiime2.core.testing')
        self.assertEqual(plugin.citations, ())

    def test_user_support_text_default(self):
        plugin = qiime2.plugin.Plugin(
            name='local-dummy-plugin',
            version='0.0.0-dev',
            website='https://github.com/qiime2/qiime2',
            package='qiime2.core.testing')
        self.assertTrue(plugin.user_support_text.startswith('Please post'))
        self.assertTrue(plugin.user_support_text.endswith(
                            'https://forum.qiime2.org'))

    def test_actions(self):
        actions = self.plugin.actions

        self.assertIsInstance(actions, types.MappingProxyType)
        self.assertEqual(actions.keys(),
                         {'merge_mappings', 'concatenate_ints', 'split_ints',
                          'most_common_viz', 'mapping_viz',
                          'identity_with_metadata',
                          'identity_with_metadata_column',
                          'identity_with_categorical_metadata_column',
                          'identity_with_numeric_metadata_column',
                          'identity_with_optional_metadata',
                          'identity_with_optional_metadata_column',
                          'params_only_method', 'no_input_method',
                          'optional_artifacts_method', 'variadic_input_method',
                          'params_only_viz', 'no_input_viz',
                          'long_description_method', 'parameter_only_pipeline',
                          'typical_pipeline', 'optional_artifact_pipeline',
                          'pointless_pipeline', 'visualizer_only_pipeline',
                          'pipelines_in_pipeline', 'failing_pipeline',
                          'docstring_order_method',
                          'constrained_input_visualization',
                          'combinatorically_mapped_method',
                          'double_bound_variable_method',
                          'bool_flag_swaps_output_method',
                          'predicates_preserved_method',
                          'deprecated_method', 'union_inputs',
                          'unioned_primitives',
                          'type_match_list_and_set',
                          })
        for action in actions.values():
            self.assertIsInstance(action, qiime2.sdk.Action)

        # Read-only dict.
        with self.assertRaises(TypeError):
            actions["i-shouldn't-do-this"] = "my-action"
        with self.assertRaises(TypeError):
            actions["merge_mappings"] = "my-action"

    def test_methods(self):
        methods = self.plugin.methods

        self.assertEqual(methods.keys(),
                         {'merge_mappings', 'concatenate_ints', 'split_ints',
                          'identity_with_metadata',
                          'identity_with_metadata_column',
                          'identity_with_categorical_metadata_column',
                          'identity_with_numeric_metadata_column',
                          'identity_with_optional_metadata',
                          'identity_with_optional_metadata_column',
                          'params_only_method', 'no_input_method',
                          'optional_artifacts_method',
                          'long_description_method',
                          'variadic_input_method', 'docstring_order_method',
                          'combinatorically_mapped_method',
                          'double_bound_variable_method',
                          'bool_flag_swaps_output_method',
                          'predicates_preserved_method',
                          'deprecated_method', 'union_inputs',
                          'unioned_primitives',
                          'type_match_list_and_set',
                          })
        for method in methods.values():
            self.assertIsInstance(method, qiime2.sdk.Method)

    def test_visualizers(self):
        visualizers = self.plugin.visualizers

        self.assertEqual(visualizers.keys(),
                         {'most_common_viz', 'mapping_viz', 'params_only_viz',
                          'no_input_viz', 'constrained_input_visualization'})
        for viz in visualizers.values():
            self.assertIsInstance(viz, qiime2.sdk.Visualizer)

    def test_pipelines(self):
        pipelines = self.plugin.pipelines

        self.assertEqual(pipelines.keys(),
                         {'parameter_only_pipeline', 'typical_pipeline',
                          'optional_artifact_pipeline', 'pointless_pipeline',
                          'visualizer_only_pipeline', 'pipelines_in_pipeline',
                          'failing_pipeline'})
        for pipeline in pipelines.values():
            self.assertIsInstance(pipeline, qiime2.sdk.Pipeline)

    # TODO test registration of directory formats.

    def test_type_fragments(self):
        types = self.plugin.type_fragments.keys()

        self.assertEqual(
            set(types),
            set(['IntSequence1', 'IntSequence2', 'IntSequence3', 'Mapping',
                 'FourInts', 'Kennel', 'Dog', 'Cat', 'SingleInt', 'C1', 'C2',
                 'C3', 'Foo', 'Bar', 'Baz']))

    def test_types(self):
        types = self.plugin.types
        # Get just the SemanticTypeRecords out of the types dictionary, then
        # get just the types out of the SemanticTypeRecord namedtuples
        types = {type_.semantic_type for type_ in types.values()}

        exp = {IntSequence1, IntSequence2, FourInts, Mapping, Kennel[Dog],
               Kennel[Cat], SingleInt}
        self.assertLessEqual(exp, types)
        self.assertNotIn(Cat, types)
        self.assertNotIn(Dog, types)
        self.assertNotIn(Kennel, types)


if __name__ == '__main__':
    unittest.main()
