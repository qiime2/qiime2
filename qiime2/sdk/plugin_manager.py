# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import os
import pkg_resources
import qiime2.core.type
from qiime2.core.util import get_view_name


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
        self._importable = set()
        self._exportable = set()
        self._canonical_formats = set()

        # These are all dependent loops, each requires the loop above it to
        # be completed.
        for entry_point in self.iter_entry_points():
            plugin = entry_point.load()
            self.plugins[plugin.name] = plugin

        for plugin in self.plugins.values():
            self._integrate_plugin(plugin)

        # TODO: Delete this loop once the helper methods for get formats are
        # approved
        for canonical_format in self._canonical_formats:
            self._importable.update(
                                self._reverse_transformers[canonical_format])
            self._exportable.update(self.transformers[canonical_format])

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
                self._ff_to_sfdf[fmt.file.format] = fmt

            # TODO: remove this when `sniff` is removed
            if hasattr(fmt, 'sniff') and hasattr(fmt, '_validate_'):
                raise RuntimeError(
                    'Format %r registered in plugin %r defines sniff and'
                    '_validate_ methods - only one is permitted.' %
                    (name, record.plugin.name)
                )

            self.formats[name] = record
        self.type_formats.extend(plugin.type_formats)

        for type_format in plugin.type_formats:
            self._canonical_formats.add(type_format.format)
            if isinstance(type_format.format,
                          qiime2.plugin.model.SingleFileDirectoryFormatBase):
                self._canonical_formats.add(type_format.format.file.format)

    def get_semantic_types(self):
        types = set()

        for plugin in self.plugins.values():
            for type_record in plugin.types:
                types.add(type_record.semantic_type)

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

        semantic_type : TypeExpression
            The semantic type is used to filter the formats associated with
            that specific semantic type

        This method will filter out the formats using the filter provided by
        the user and the semantic type. The return is a set of filtered
        formats.
        """

        if semantic_type is not None and semantic_type \
                not in self.get_semantic_types():
            raise ValueError("The semantic type provided: %s is not "
                             "valid.", (semantic_type))

        if filter is not None and filter not in GetFormatFilters:
            raise ValueError("The format filter provided: %s is not "
                             "valid.", (filter))

        if semantic_type is None:
            formats = self.formats

        else:
            formats = {}
            # loop through the canonical formats
            for type_format in self.type_formats:
                # if the semantic_type param is in the type exp in for a
                # canonical format, then loop through the format, record in
                # formats
                if semantic_type <= type_format.type_expression:
                    for single_format, record in self.formats:
                        # if the format in formats is the format in canonical,
                        # then set the format as the key and the record as the
                        # value (record is a formatrecord type)
                        # the return is a formatrecord dict
                        if self.formats[single_format].format is \
                                type_format.format:
                            formats[type_format.format] = record

        # These variables contain format -> transrecord dict
        importable_formats = self._get_importable_formats(self.formats)
        exportable_formats = self._get_exportable_formats(self.formats)

        # The filtering, using set comprehension, is comparing formatrecs to
        # transrecords if the filter is none
        if filter is None:
            formats = {**formats, **importable_formats, **exportable_formats}
        elif filter is GetFormatFilters.IMPORTABLE:
            formats = importable_formats
        elif filter is GetFormatFilters.EXPORTABLE:
            formats = exportable_formats
        elif filter is GetFormatFilters.IMPORTABLE | \
                GetFormatFilters.EXPORTABLE:
            formats = {**importable_formats, **exportable_formats}

        # The return is either a formatrecord dict or transrecord dict
        return formats

    # This method takes in an iterable of formats to perform looping in helper
    def _get_exportable_formats(self, formats):
        """
        _get_exportable_formats(self, formats)

        formats : DirectoryFormat
            We are finding all formats that can be directly transformed to from
            this format

        This method creates a set utilizing the transformers dictionary and
        the canonical format passed to it as a key. This allows
        for access to formats the canonical format transforms into.
        """
        exportable_formats = {}
        temp_list = {}
        SFDF = qiime2.plugin.model.SingleFileDirectoryFormatBase

        for format_, record in formats.items():
            # Need to get the class of the format record to check against the
            # parent class
            if issubclass(record.format, SFDF):
                format_name = get_view_name(record.format.file.format)
                # TODO: Need to build format record for FileFormat (Possible)
                temp_list[format_name] = record.format.file.format

                # temp_list[format_.format] = format_.format.file.format

            temp_list[format_] = formats[format_]
            # temp_list.append(formats[format_name])

        for format_item in temp_list:

            if format_item in self.transformers:
                # exportable_formats.extend(self.transformers[format_item])
                exportable_formats[format_item] = \
                    self.transformers[format_item]

                for format_dict_item in self.transformers[format_item]:

                    if issubclass(format_dict_item, SFDF):
                        # exportable_formats.append(format_dict_item.file.format)
                        format_name = get_view_name(
                            format_dict_item.file.format)
                        temp_list[format_name] = formats[format_name]

                    if format_dict_item in self._ff_to_sfdf:
                        exportable_formats[format_dict_item] = \
                            self._ff_to_sfdf[format_dict_item]

        exportable_formats.update(temp_list)

        return exportable_formats

    # This method takes in an iterable of formats to perform looping in helper
    def _get_importable_formats(self, formats):
        """
        _get_importable_formats(self, formats)

        formats : DirectoryFormat
            We are finding all formats that can be directly transformed into
            this format

        This method creates a set utilizing the reverse transformers
        dictionary and the canonical format passed to it as a key. This allows
        for access to formats that transform into the canonical format.
        """
        importable_formats = {}
        SFDF = qiime2.plugin.model.SingleFileDirectoryFormatBase

        temp_list = {}
        SFDF = qiime2.plugin.model.SingleFileDirectoryFormatBase

        for format_, record in formats.items():
            # Need to get the class of the format record to check against the
            # parent class
            if issubclass(record.format, SFDF):
                format_name = get_view_name(record.format.file.format)
                # TODO: Need to build format record for FileFormat (Possible)
                temp_list[format_name] = record.format.file.format

                # temp_list[format_.format] = format_.format.file.format

            temp_list[format_] = formats[format_]
            # temp_list.append(formats[format_name])

        for format_item in temp_list:

            if format_item in self._reverse_transformers:
                # exportable_formats.extend(self.transformers[format_item])
                importable_formats[format_item] = \
                    self.transformers[format_item]

                for format_dict_item in \
                        self._reverse_transformers[format_item]:

                    if issubclass(format_dict_item, SFDF):
                        # exportable_formats.append(format_dict_item.file.format)
                        format_name = get_view_name(
                            format_dict_item.file.format)
                        temp_list[format_name] = formats[format_name]

                    if format_dict_item in self._ff_to_sfdf:
                        importable_formats[format_dict_item] = \
                            self._ff_to_sfdf[format_dict_item]

        importable_formats.update(temp_list)

        '''
        for format_ in formats:

            if issubclass(format_, SFDF):
                importable_formats.append(format_.file.format)

            importable_formats.extend(format_)

        for format_item in importable_formats:

            if format_item in self.transformers:
                importable_formats.extend(self.transformers[format_item])

                for format_dict_item in self.transformers[format_item]:

                    if issubclass(format_dict_item, SFDF):
                        importable_formats.append(format_dict_item)

                    if format_dict_item in self._ff_to_sfdf:
                        importable_formats.append(format_dict_item)
        '''
        return importable_formats

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
