# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
import unittest
import unittest.mock

import qiime.core.testing
from qiime.plugin import Plugin
from qiime.sdk import PluginManager


class TestPluginManager(unittest.TestCase):
    def test_constructor(self):
        class MockEntryPoint(pkg_resources.EntryPoint):
            def __init__(self, plugin):
                self.plugin = plugin

            def load(self):
                return self.plugin

        plugin1 = Plugin(
            name='dummy-plugin',
            version='0.4.2',
            website='www.dummy-plugin-hub.com',
            package='dummy_plugin'
        )
        plugin2 = Plugin(
            name='other-dummy-plugin',
            version='0.4.2',
            website='www.dummy-plugin-hub.com',
            package='dummy_plugin'
        )

        def mock_iter_entry_points(group):
            return [MockEntryPoint(plugin1), MockEntryPoint(plugin2)]

        with unittest.mock.patch.object(pkg_resources, 'iter_entry_points',
                                        mock_iter_entry_points):
            plugin_manager = PluginManager()

        self.assertEqual(
            plugin_manager.plugins,
            {'dummy-plugin': plugin1, 'other-dummy-plugin': plugin2,
             'test-plugin': qiime.core.testing.plugin})

        expected_archive_format = qiime.plugin.ArchiveFormat(
            'example-archive-format', 1, lambda e: e)
        self.assertEqual(
            plugin_manager.archive_formats,
            {('example-archive-format', 1):
             ('test-plugin', expected_archive_format)})


if __name__ == '__main__':
    unittest.main()
