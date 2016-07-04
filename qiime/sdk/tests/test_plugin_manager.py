# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import unittest

import qiime.plugin
import qiime.sdk
from qiime.core.data_layout import DataLayout

from qiime.core.testing.data_layout import (int_sequence_to_list,
                                            int_sequence_to_counter,
                                            list_to_int_sequence,
                                            mapping_to_dict,
                                            dict_to_mapping)
from qiime.core.testing.type import (IntSequence1, IntSequence2, Mapping,
                                     FourInts)
from qiime.core.testing.util import get_dummy_plugin


# These tests are written to pass regardless of other plugins that may be
# installed. This is accomplished by only testing pieces we know are defined by
# the dummy plugin.
class TestPluginManager(unittest.TestCase):
    def setUp(self):
        # Ignore the returned dummy plugin object, just run this to verify the
        # plugin exists as the tests rely on it being loaded. PluginManager is
        # a singleton so there's no issue creating it again.
        get_dummy_plugin()
        self.pm = qiime.sdk.PluginManager()

    def test_plugins(self):
        plugins = self.pm.plugins

        self.assertIsInstance(plugins['dummy-plugin'], qiime.plugin.Plugin)

    def test_data_layouts(self):
        data_layouts = self.pm.data_layouts

        self.assertIn(('int-sequence', 1), data_layouts)
        self.assertIn(('mapping', 1), data_layouts)
        self.assertIn(('four-ints', 1), data_layouts)

        plugin_name, data_layout = data_layouts[('int-sequence', 1)]

        self.assertEqual(plugin_name, 'dummy-plugin')

        self.assertIsInstance(data_layout, DataLayout)
        self.assertEqual(data_layout.name, 'int-sequence')
        self.assertEqual(data_layout.version, 1)

        self.assertEqual(data_layout.readers,
                         {list: int_sequence_to_list,
                          collections.Counter: int_sequence_to_counter})
        self.assertEqual(data_layout.writers, {list: list_to_int_sequence})

        with self.assertRaises(NotImplementedError):
            data_layout.validate('some-data-directory')

    def test_semantic_types(self):
        types = self.pm.semantic_types

        self.assertEqual(types['IntSequence1'], ('dummy-plugin', IntSequence1))
        self.assertEqual(types['IntSequence2'], ('dummy-plugin', IntSequence2))
        self.assertEqual(types['Mapping'], ('dummy-plugin', Mapping))
        self.assertEqual(types['FourInts'], ('dummy-plugin', FourInts))

    def test_semantic_type_to_data_layouts(self):
        type_to_layouts = self.pm.semantic_type_to_data_layouts

        self.assertIsInstance(type_to_layouts[IntSequence1], DataLayout)
        self.assertIsInstance(type_to_layouts[IntSequence2], DataLayout)
        self.assertIsInstance(type_to_layouts[Mapping], DataLayout)
        self.assertIsInstance(type_to_layouts[FourInts], DataLayout)

        data_layout = type_to_layouts[Mapping]

        self.assertEqual(data_layout.name, 'mapping')
        self.assertEqual(data_layout.version, 1)

        self.assertEqual(data_layout.readers, {dict: mapping_to_dict})
        self.assertEqual(data_layout.writers, {dict: dict_to_mapping})

        with self.assertRaises(NotImplementedError):
            data_layout.validate('some-data-directory')


if __name__ == '__main__':
    unittest.main()
