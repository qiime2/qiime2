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
import qiime2.core.type
from qiime2.core.format import FormatBase
from qiime2.plugin.model import SingleFileDirectoryFormatBase
from qiime2.sdk.util import parse_type

import enum


class GetFormatFilters(enum.Flag):
    EXPORTABLE = enum.auto()
    IMPORTABLE = enum.auto()


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
    def __new__(cls):
        if cls.__instance is None:
            self = super().__new__(cls)
            self._init()
            cls.__instance = self
        return cls.__instance

    def _init(self):
        self.plugins = {}
        self.type_fragments = {}
        self.transformers = collections.defaultdict(dict)
        self._reverse_transformers = collections.defaultdict(dict)
        self.formats = {}
        self.views = {}
        self.type_formats = []
        self._ff_to_sfdf = {}

        # These are all dependent loops, each requires the loop above it to
        # be completed.
        for entry_point in self.iter_entry_points():
            plugin = entry_point.load()
            self.plugins[plugin.name] = plugin

        for plugin in self.plugins.values():
            self._integrate_plugin(plugin)

    def _integrate_plugin(self, plugin):
        for type_name, type_record in plugin.type_fragments.items():
            if type_name in self.type_fragments:
                conflicting_type_record = \
                    self.type_fragments[type_name]
                raise ValueError("Duplicate semantic type (%r) defined in"
                                 " plugins: %r and %r"
                                 % (type_name, type_record.plugin.name,
                                    conflicting_type_record.plugin.name))

            self.type_fragments[type_name] = type_record

        for (input, output), transformer_record in plugin.transformers.items():
            if output in self.transformers[input]:
                raise ValueError("Transformer from %r to %r already exists."
                                 % transformer_record)
            self.transformers[input][output] = transformer_record
            self._reverse_transformers[output][input] = transformer_record

        for name, record in plugin.views.items():
            if name in self.views:
                raise NameError(
                    "Duplicate view registration (%r) defined in plugins: %r"
                    " and %r" %
                    (name, record.plugin.name, self.formats[name].plugin.name)
                )
            self.views[name] = record

        for name, record in plugin.formats.items():
            fmt = record.format

            if issubclass(
                    fmt, qiime2.plugin.model.SingleFileDirectoryFormatBase):
                if fmt.file.format in self._ff_to_sfdf.keys():
                    self._ff_to_sfdf[fmt.file.format].add(fmt)
                else:
                    self._ff_to_sfdf[fmt.file.format] = {fmt}

            # TODO: remove this when `sniff` is removed
            if hasattr(fmt, 'sniff') and hasattr(fmt, '_validate_'):
                raise RuntimeError(
                    'Format %r registered in plugin %r defines sniff and'
                    '_validate_ methods - only one is permitted.' %
                    (name, record.plugin.name)
                )

            self.formats[name] = record
        self.type_formats.extend(plugin.type_formats)

    def get_semantic_types(self):
        types = set()

        for plugin in self.plugins.values():
            for type_record in plugin.types:
                types.add(type_record)

        return types

    # TODO: Should plugin loading be transactional? i.e. if there's
    # something wrong, the entire plugin fails to load any piece, like a
    # databases rollback/commit

    def get_formats(self, *, filter=None, semantic_type=None):
        """
        get_formats(self, *, filter=None, semantic_type=None)

        filter : enum
            filter is an enum integer that will be used to determine user
            input to output specified formats

        semantic_type : TypeExpression | String
            The semantic type is used to filter the formats associated with
            that specific semantic type

        This method will filter out the formats using the filter provided by
        the user and the semantic type. The return is a set of filtered
        formats.
        """
        if filter is not None and filter not in GetFormatFilters:
            raise ValueError("The format filter provided: %s is not "
                             "valid.", (filter))

        if semantic_type is None:
            formats = set(f.format for f in self.type_formats)

        else:
            formats = set()

            if isinstance(semantic_type, str):
                semantic_type = parse_type(semantic_type, "semantic")

            for type_format in self.type_formats:
                if semantic_type <= type_format.type_expression:
                    formats.add(type_format.format)
                    break

            if not formats:
                raise ValueError("No formats associated with the type "
                                 f"{semantic_type}.")

        transformable_formats = set(formats)

        if filter is None or GetFormatFilters.IMPORTABLE in filter:
            transformable_formats.update(
                self._get_formats_helper(formats, self._reverse_transformers))

        if filter is None or GetFormatFilters.EXPORTABLE in filter:
            transformable_formats.update(
                self._get_formats_helper(formats, self.transformers))

        result_formats = {}
        for format_ in transformable_formats:
            format_ = format_.__name__
            result_formats[format_] = self.formats[format_]

        return result_formats

    def _get_formats_helper(self, formats, transformer_dict):
        """
        _get_formats_helper(self, formats, transformer_dict)

        formats : set(DirectoryFormat)
            We are finding all formats that can be directly transformed to and
            from formats in this set

        tranformer_dict : Class: Dict{ Class: TransformerRecord }
            The dictionary of transformers allows the method to get formats
            that are transformable from the given format

        This method creates a set utilizing the transformers dictionary and
        the formats set to get related formats for a specific format.
        """
        query_set = set(formats)

        for format_ in formats:
            if issubclass(format_, SingleFileDirectoryFormatBase):
                if format_.file.format.__name__ in self.formats:
                    query_set.add(format_.file.format)

        result_formats = set(query_set)

        for format_ in query_set:
            for transformed_format in transformer_dict[format_]:
                if issubclass(transformed_format, FormatBase):
                    result_formats.add(transformed_format)

                    if issubclass(transformed_format,
                                  SingleFileDirectoryFormatBase):
                        result_formats.add(transformed_format.file.format)

                    if transformed_format in self._ff_to_sfdf:
                        result_formats.update(
                            self._ff_to_sfdf[transformed_format])

        return result_formats

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
