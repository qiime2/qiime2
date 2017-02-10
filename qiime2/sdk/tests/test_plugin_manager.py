# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import qiime2.plugin
import qiime2.sdk
from qiime2.plugin.plugin import SemanticTypeRecord

from qiime2.core.testing.type import (IntSequence1, IntSequence2, Mapping,
                                      FourInts, Kennel, Dog, Cat)
from qiime2.core.testing.util import get_dummy_plugin


class TestPluginManager(unittest.TestCase):
    def setUp(self):
        self.plugin = get_dummy_plugin()
        # PluginManager is a singleton so there's no issue creating it again.
        self.pm = qiime2.sdk.PluginManager()

    def test_plugins(self):
        plugins = self.pm.plugins

        exp = {'dummy-plugin': self.plugin}
        self.assertEqual(plugins, exp)

    def test_semantic_types(self):
        types = self.pm.semantic_types

        exp = {
            'IntSequence1': SemanticTypeRecord(semantic_type=IntSequence1,
                                               plugin=self.plugin),
            'IntSequence2': SemanticTypeRecord(semantic_type=IntSequence2,
                                               plugin=self.plugin),
            'Mapping': SemanticTypeRecord(semantic_type=Mapping,
                                          plugin=self.plugin),
            'FourInts': SemanticTypeRecord(semantic_type=FourInts,
                                           plugin=self.plugin),
            'Kennel': SemanticTypeRecord(semantic_type=Kennel,
                                         plugin=self.plugin),
            'Dog': SemanticTypeRecord(semantic_type=Dog,
                                      plugin=self.plugin),
            'Cat': SemanticTypeRecord(semantic_type=Cat,
                                      plugin=self.plugin),
        }

        self.assertEqual(types, exp)

    def test_importable_types(self):
        types = self.pm.importable_types

        exp = {IntSequence1, IntSequence2, FourInts, Mapping, Kennel[Dog],
               Kennel[Cat]}
        self.assertEqual(types, exp)

    # TODO: add tests for type/directory/transformer registrations


if __name__ == '__main__':
    unittest.main()
