# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import os
import pkg_resources

import qiime.core.type


class PluginManager:
    __instance = None

    # This class is a singleton as it is slow to create, represents the
    # state of a qiime installation, and is needed *everywhere*
    def __new__(cls):
        if cls.__instance is None:
            self = super().__new__(cls)
            self._init()
            cls.__instance = self
        return cls.__instance

    def _init(self):
        self.plugins = {}
        self.semantic_types = {}
        self.transformers = collections.defaultdict(dict)
        self.type_formats = []

        plugins = []
        for entry_point in pkg_resources.iter_entry_points(
                group='qiime.plugins'):
            if entry_point.name != 'dummy-plugin' or 'QIIMETEST' in os.environ:
                plugins.append(entry_point.load())

        # These are all dependent loops, each requires the loop above it to
        # be completed.
        for plugin in plugins:
            self.plugins[plugin.name] = plugin

        for plugin in plugins:
            self._integrate_plugin(plugin)

    def _integrate_plugin(self, plugin):
        for type_name, type_record in plugin.types.items():
            if type_name in self.semantic_types:
                conflicting_type_record = self.semantic_types[type_name]
                raise ValueError("Duplicate semantic type (%r) defined in"
                                 " plugins: %r and %r"
                                 % (type_name, type_record.plugin.name,
                                    conflicting_type_record.plugin.name))

            self.semantic_types[type_name] = type_record

        for (input, output), transformer_record in plugin.transformers.items():
            if output in self.transformers[input]:
                raise ValueError("Transformer from %r to %r already exists."
                                 % transformer_record)
            self.transformers[input][output] = transformer_record

        self.type_formats.extend(plugin.type_formats)

    # TODO: Should plugin loading be transactional? i.e. if there's
    # something wrong, the entire plugin fails to load any piece, like a
    # databases rollback/commit

    def get_directory_format(self, semantic_type):
        if not qiime.core.type.is_semantic_type(semantic_type):
            raise TypeError(
                "Must provide a semantic type via `type`, not %r" % type)

        dir_fmt = None
        for type_format_record in self.type_formats:
            if semantic_type <= type_format_record.type_expression:
                dir_fmt = type_format_record.format
                break

        if dir_fmt is None:
            raise TypeError(
                "Semantic type %r does not have a compatible directory format."
                % semantic_type)

        return dir_fmt
