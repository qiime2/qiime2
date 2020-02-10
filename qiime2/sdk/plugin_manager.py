# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import os
import pkg_resources
import qiime2.core.transform as transform
import qiime2.core.type


class PluginManager:
    entry_point_group = 'qiime2.plugins'
    __instance = None

    @classmethod
    def iter_entry_points(cls):
        """Yield QIIME 2 plugin entry points.

        If the QIIMETEST environment variable is set, only the framework
        testing plugin entry point (`dummy-plugin`) will be yielded. Otherwise,
        all available plugin entry points (excluding `dummy-plugin`) will be
        yielded.

        """
        for entry_point in pkg_resources.iter_entry_points(
                group=cls.entry_point_group):
            if 'QIIMETEST' in os.environ:
                if entry_point.name == 'dummy-plugin':
                    yield entry_point
            else:
                if entry_point.name != 'dummy-plugin':
                    yield entry_point

    # This class is a singleton as it is slow to create, represents the
    # state of a qiime2 installation, and is needed *everywhere*
    def __new__(cls, add_plugins=True):
        if cls.__instance is None:
            self = super().__new__(cls)
            self._init(add_plugins=add_plugins)
            cls.__instance = self
        else:
            if add_plugins is False:
                raise ValueError(
                    'PluginManager singleton already exists, cannot change '
                    'default value for `add_plugins`.')
        return cls.__instance

    def _init(self, add_plugins):
        self.plugins = {}
        self._plugin_by_id = {}
        self.semantic_types = {}
        self.transformers = collections.defaultdict(dict)
        self.formats = {}
        self.views = {}
        self.type_formats = []

        if add_plugins:
            # These are all dependent loops, each requires the loop above it to
            # be completed.
            for entry_point in self.iter_entry_points():
                project_name = entry_point.dist.project_name
                package = entry_point.module_name.split('.')[0]
                plugin = entry_point.load()

                self.add_plugin(plugin, package, project_name)

    def add_plugin(self, plugin, package=None, project_name=None):
        self.plugins[plugin.name] = plugin
        self._plugin_by_id[plugin.id] = plugin
        if plugin.package is None:
            plugin.package = package
        if plugin.project_name is None:
            plugin.project_name = project_name

        # validate _after_ applying arguments
        if plugin.package is None:
            raise ValueError(
                'No value specified for package - must provide a value for '
                '`package` or set `plugin.package`.')
        if plugin.project_name is None:
            raise ValueError(
                'No value specified for project_name - must proved a value '
                'for `project_name` or set `plugin.project_name`.')

        self._integrate_plugin(plugin)
        plugin.freeze()

    def get_plugin(self, *, id=None, name=None):
        if id is None and name is None:
            raise ValueError("No plugin requested.")
        elif id is not None:
            try:
                return self._plugin_by_id[id]
            except KeyError:
                raise KeyError('No plugin currently registered '
                               'with id: "%s".' % (id,))
        else:
            try:
                return self.plugins[name]
            except KeyError:
                raise KeyError('No plugin currently registered '
                               'with name: "%s".' % (name,))

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

        for name, record in plugin.views.items():
            if name in self.views:
                raise NameError(
                    "Duplicate view registration (%r) defined in plugins: %r"
                    " and %r" %
                    (name, record.plugin.name, self.formats[name].plugin.name)
                )
            self.views[name] = record

        for name, record in plugin.formats.items():
            # TODO: remove this when `sniff` is removed
            fmt = record.format
            if hasattr(fmt, 'sniff') and hasattr(fmt, '_validate_'):
                raise RuntimeError(
                    'Format %r registered in plugin %r defines sniff and'
                    '_validate_ methods - only one is permitted.' %
                    (name, record.plugin.name)
                )

            self.formats[name] = record
        self.type_formats.extend(plugin.type_formats)

    # TODO: Should plugin loading be transactional? i.e. if there's
    # something wrong, the entire plugin fails to load any piece, like a
    # databases rollback/commit

    @property
    def importable_formats(self):
        """Return formats that are importable.

        A format is importable in a QIIME 2 deployment if it can be transformed
        into at least one of the canonical semantic type formats.

        """
        importable_formats = {}
        for name, record in self.formats.items():
            from_type = transform.ModelType.from_view_type(
                record.format)
            for type_format in self.type_formats:
                to_type = transform.ModelType.from_view_type(
                    type_format.format)
                if from_type.has_transformation(to_type):
                    importable_formats[name] = record
                    break
        return importable_formats

    @property
    def importable_types(self):
        """Return set of concrete semantic types that are importable.

        A concrete semantic type is importable if it has an associated
        directory format.

        """
        importable_types = set()
        for type_format in self.type_formats:
            for type in type_format.type_expression:
                importable_types.add(type)
        return importable_types

    def get_directory_format(self, semantic_type):
        if not qiime2.core.type.is_semantic_type(semantic_type):
            raise TypeError(
                "Must provide a semantic type via `semantic_type`, not %r" %
                semantic_type)

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
