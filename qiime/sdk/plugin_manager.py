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
        self.semantic_type_to_archive_formats = {}

        plugins = [e.load() for e in pkg_resources.iter_entry_points(
            group='qiime.plugin')]
        if '__QIIME2_TEST_MODE__' in os.environ:
            import qiime.core.testing
            plugins.append(qiime.core.testing.plugin)

        # These are all dependent loops, each requires the loop above it to
        # be completed.
        for plugin in plugins:
            self.plugins[plugin.name] = plugin

        for plugin in plugins:
            self._integrate_plugin(plugin)

        for plugin in plugins:
            self._finalize_plugin(plugin)

    def _integrate_plugin(self, plugin):
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

    def _finalize_plugin(self, plugin):
        for semantic_type, archive_format_id in \
                plugin.type_to_archive_formats.items():
            if archive_format_id not in self.archive_formats:
                raise ValueError("Archive format %r does not exist, cannot"
                                 " register semantic type (%r) to it."
                                 % (archive_format_id, semantic_type))
            self.semantic_type_to_archive_formats[semantic_type] = \
                self.archive_formats[archive_format_id][1]

        for reader_registration, reader in \
                plugin.archive_format_readers.items():
            name, version, view_type = reader_registration
            _, archive_format = self.archive_formats[(name, version)]
            archive_format.reader(view_type)(reader)

        for writer_registration, writer in \
                plugin.archive_format_writers.items():
            name, version, view_type = writer_registration
            _, archive_format = self.archive_formats[(name, version)]
            archive_format.writer(view_type)(writer)

    # TODO: Should plugin loading be transactional? i.e. if there's
    # something wrong, the entire plugin fails to load any piece, like a
    # databases rollback/commit
