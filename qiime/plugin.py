# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources

import qiime.sdk
import qiime.core.data_layout
import qiime.core.type.grammar as grammar
from qiime.core.type import (SemanticType, is_semantic_type, Int, Str, Float,
                             Color, Metadata, MetadataCategory, Properties,
                             Range, Choices, Arguments)

__all__ = ['Plugin', 'Int', 'Str', 'Float', 'Color', 'Metadata',
           'MetadataCategory', 'SemanticType', 'Properties', 'Range',
           'Choices', 'Arguments']


class Plugin:
    def __init__(self, name, version, website, package):
        self.name = name
        self.version = version
        self.website = website
        self.package = package

        self.methods = PluginMethods(self.package)
        self.visualizers = PluginVisualizers()

        self.data_layouts = {}
        self.types = {}
        self.type_to_data_layouts = {}
        self.data_layout_readers = {}
        self.data_layout_writers = {}

    def register_data_layout_reader(self, name, version, view_type, reader):
        self.data_layout_readers[(name, version, view_type)] = reader

    def register_data_layout_writer(self, name, version, view_type, writer):
        self.data_layout_writers[(name, version, view_type)] = writer

    def register_data_layout(self, name, version, validator):
        self.data_layouts[(name, version)] = qiime.core.data_layout.DataLayout(
            name, version, validator)

    def register_semantic_type(self, semantic_type):
        if not is_semantic_type(semantic_type):
            raise TypeError("%r is not a semantic type." % semantic_type)

        if not (isinstance(semantic_type, grammar.CompositeType) or
                (semantic_type.is_concrete() and not semantic_type.fields)):
            raise ValueError("%r is not a semantic type symbol."
                             % semantic_type)

        if semantic_type.name in self.types:
            raise ValueError("Duplicate semantic type symbol %r."
                             % semantic_type)

        self.types[semantic_type.name] = semantic_type

    def register_type_to_data_layout(self, semantic_type, name, version):
        if not is_semantic_type(semantic_type):
            raise TypeError("%r is not a semantic type." % semantic_type)
        if not isinstance(semantic_type, grammar.TypeExpression):
            raise ValueError("%r is not a semantic type expression."
                             % semantic_type)

        self.type_to_data_layouts[semantic_type] = (name, version)

    # TODO: should register_data_layout, register_semantic_type, and
    # register_type_to_data_layout be namespaced to indicate that they
    # shouldn't typically be used by plugin devs?


# TODO refactor these classes, they are basically the same
class PluginMethods(dict):
    def __init__(self, package):
        self._package = package
        super().__init__()

    def register_function(self, function, inputs, parameters, outputs, name,
                          description):
        method = qiime.sdk.Method.from_function(function, inputs, parameters,
                                                outputs, name, description)
        self[method.id] = method

    def register_markdown(self, markdown_filepath):
        markdown_filepath = pkg_resources.resource_filename(
            self._package, markdown_filepath)
        method = qiime.sdk.Method.from_markdown(markdown_filepath)
        self[method.id] = method


class PluginVisualizers(dict):
    def register_function(self, function, inputs, parameters, name,
                          description):
        visualizer = qiime.sdk.Visualizer.from_function(
            function, inputs, parameters, name, description)
        self[visualizer.id] = visualizer
