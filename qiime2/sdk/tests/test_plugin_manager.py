# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import qiime2.plugin
import qiime2.sdk

from qiime2.core.testing.type import (IntSequence1, IntSequence2, Mapping,
                                      FourInts)
from qiime2.core.testing.util import get_dummy_plugin


# These tests are written to pass regardless of other plugins that may be
# installed. This is accomplished by only testing pieces we know are defined by
# the dummy plugin.
class TestPluginManager(unittest.TestCase):
    def setUp(self):
        self.plugin = get_dummy_plugin()
        # PluginManager is a singleton so there's no issue creating it again.
        self.pm = qiime2.sdk.PluginManager()

    def test_plugins(self):
        plugins = self.pm.plugins

        self.assertIsInstance(plugins['dummy-plugin'], qiime2.plugin.Plugin)

    def test_semantic_types(self):
        types = self.pm.semantic_types

        self.assertEqual(types['IntSequence1'].semantic_type, IntSequence1)
        self.assertEqual(types['IntSequence1'].plugin, self.plugin)

        self.assertEqual(types['IntSequence2'].semantic_type, IntSequence2)
        self.assertEqual(types['IntSequence2'].plugin, self.plugin)

        self.assertEqual(types['Mapping'].semantic_type, Mapping)
        self.assertEqual(types['Mapping'].plugin, self.plugin)

        self.assertEqual(types['FourInts'].semantic_type, FourInts)
        self.assertEqual(types['FourInts'].plugin, self.plugin)

    # TODO: add tests for type/directory/transformer registrations


if __name__ == '__main__':
    unittest.main()
