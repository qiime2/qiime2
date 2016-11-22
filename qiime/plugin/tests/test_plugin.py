# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import types
import unittest

import qiime.plugin
import qiime.sdk

from qiime.core.testing.util import get_dummy_plugin


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
        self.assertEqual(self.plugin.package, 'qiime.core.testing')

    def test_citation_text(self):
        self.assertEqual(self.plugin.citation_text, 'No relevant citation.')

    def test_user_support_text(self):
        self.assertEqual(self.plugin.user_support_text,
                         'For help, see https://qiime2.org')

    def test_citation_text_default(self):
        plugin = qiime.plugin.Plugin(
            name='local-dummy-plugin',
            version='0.0.0-dev',
            website='https://github.com/qiime2/qiime2',
            package='qiime.core.testing')
        self.assertTrue(plugin.citation_text.startswith('No citation'))
        self.assertTrue(plugin.citation_text.endswith(
                            plugin.website))

    def test_user_support_text_default(self):
        plugin = qiime.plugin.Plugin(
            name='local-dummy-plugin',
            version='0.0.0-dev',
            website='https://github.com/qiime2/qiime2',
            package='qiime.core.testing')
        self.assertTrue(plugin.user_support_text.startswith('No user'))
        self.assertTrue(plugin.user_support_text.endswith(
                            plugin.website))

    def test_actions(self):
        actions = self.plugin.actions

        self.assertIsInstance(actions, types.MappingProxyType)
        self.assertEqual(actions.keys(),
                         {'merge_mappings', 'concatenate_ints', 'split_ints',
                          'most_common_viz', 'mapping_viz',
                          'identity_with_metadata',
                          'identity_with_metadata_category'})
        for action in actions.values():
            self.assertIsInstance(action, qiime.sdk.Action)

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
                          'identity_with_metadata_category'})
        for method in methods.values():
            self.assertIsInstance(method, qiime.sdk.Action)

    def test_visualizers(self):
        visualizers = self.plugin.visualizers

        self.assertEqual(visualizers.keys(),
                         {'most_common_viz', 'mapping_viz'})
        for viz in visualizers.values():
            self.assertIsInstance(viz, qiime.sdk.Action)

    # TODO test registration of directory formats.

    def test_types(self):
        types = self.plugin.types.keys()

        self.assertEqual(
            set(types),
            set(['IntSequence1', 'IntSequence2', 'Mapping', 'FourInts']))


if __name__ == '__main__':
    unittest.main()
