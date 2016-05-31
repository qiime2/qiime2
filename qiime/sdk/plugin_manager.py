# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources


class PluginManager:
    def __init__(self):
        self.plugins = {}
        self.archive_formats = {}
        self.semantic_types = {}

        for entry_point in pkg_resources.iter_entry_points(
                group='qiime.plugin'):
            plugin = entry_point.load()
            self._load_plugin(plugin)

        if '__QIIME2_TEST_MODE__' in os.environ:
            import qiime.core.testing
            self._load_plugin(qiime.core.testing.plugin)

    def _load_plugin(self, plugin):
        self.plugins[plugin.name] = plugin

        for id_, archive_format in plugin.archive_formats.items():
            if id_ in self.archive_formats:
                conflicting_plugin, _ = self.archive_formats[id_]
                raise ValueError(
                    "Duplicate archive format defined in plugins "
                    "%r and %r: %r, %r" % (conflicting_plugin, plugin.name,
                                           id_[0], id_[1]))
            # TODO rethink the structure for mapping archive format to plugin
            # (also applies to type system and workflows).
            self.archive_formats[id_] = (plugin.name, archive_format)

        for type_name, type_expr in plugin.types.items():
            if type_name in self.semantic_types:
                conflicting_plugin, _ = self.semantic_types[type_name]
                raise ValueError("Duplicate semantic type (%r) defined in"
                                 " plugins: %r and %r"
                                 % (type_expr, plugin.name,
                                    conflicting_plugin.name))

            self.semantic_types[type_name] = (plugin.name, type_expr)

        # TODO: Should plugin loading be transactional? i.e. if there's
        # something wrong, the entire plugin fails to load any piece, like a
        # databases rollback/commit
