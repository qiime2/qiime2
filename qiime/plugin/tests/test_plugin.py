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

from qiime.core.testing.type import (IntSequence1, IntSequence2, Mapping,
                                     FourInts)
from qiime.core.testing.data_layout import (IntSequenceFormat,
                                            int_sequence_to_list,
                                            int_sequence_to_counter,
                                            list_to_int_sequence,
                                            MappingFormat,
                                            mapping_to_dict,
                                            dict_to_mapping,
                                            SingleIntFormat,
                                            four_ints_to_list,
                                            list_to_four_ints)
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

    def test_methods(self):
        methods = self.plugin.methods

        self.assertEqual(methods.keys(),
                         {'merge_mappings', 'concatenate_ints',
                          'concatenate_ints_markdown', 'split_ints',
                          'split_ints_markdown'})
        for method in methods.values():
            self.assertIsInstance(method, qiime.sdk.Method)

    def test_visualizers(self):
        visualizers = self.plugin.visualizers

        self.assertEqual(visualizers.keys(),
                         {'most_common_viz', 'mapping_viz'})
        for viz in visualizers.values():
            self.assertIsInstance(viz, qiime.sdk.Visualizer)

    def test_data_layouts(self):
        data_layouts = self.plugin.data_layouts

        self.assertEqual(
            data_layouts.keys(),
            {('int-sequence', 1), ('mapping', 1), ('four-ints', 1)})
        for dl in data_layouts.values():
            self.assertIsInstance(dl, qiime.plugin.DataLayout)

        self.assertEqual(
            data_layouts[('int-sequence', 1)].files,
            {'ints.txt': IntSequenceFormat})

        self.assertEqual(
            data_layouts[('mapping', 1)].files,
            {'mapping.tsv': MappingFormat})

        self.assertEqual(data_layouts[('four-ints', 1)].files, {
            'file1.txt': SingleIntFormat,
            'file2.txt': SingleIntFormat,
            'nested/file3.txt': SingleIntFormat,
            'nested/file4.txt': SingleIntFormat
        })

    def test_data_layouts_finalized(self):
        for data_layout in self.plugin.data_layouts.values():
            with self.assertRaisesRegex(RuntimeError, "DataLayout.*finalized"):
                data_layout.register_file('new-file.txt', SingleIntFormat)

    def test_types(self):
        types = self.plugin.types

        self.assertEqual(
            types,
            {'IntSequence1': IntSequence1, 'IntSequence2': IntSequence2,
             'Mapping': Mapping, 'FourInts': FourInts})

    def test_type_to_data_layouts(self):
        type_to_data_layouts = self.plugin.type_to_data_layouts

        self.assertEqual(type_to_data_layouts,
                         {IntSequence1: ('int-sequence', 1),
                          IntSequence2: ('int-sequence', 1),
                          Mapping: ('mapping', 1),
                          FourInts: ('four-ints', 1)})

    def test_data_layout_readers(self):
        data_layout_readers = self.plugin.data_layout_readers

        self.assertEqual(
            data_layout_readers,
            {('int-sequence', 1, list): int_sequence_to_list,
             ('int-sequence', 1, collections.Counter): int_sequence_to_counter,
             ('mapping', 1, dict): mapping_to_dict,
             ('four-ints', 1, list): four_ints_to_list})

    def test_data_layout_writers(self):
        data_layout_writers = self.plugin.data_layout_writers

        self.assertEqual(
            data_layout_writers,
            {('int-sequence', 1, list): list_to_int_sequence,
             ('mapping', 1, dict): dict_to_mapping,
             ('four-ints', 1, list): list_to_four_ints})


if __name__ == '__main__':
    unittest.main()
