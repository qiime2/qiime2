# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources

import qiime.sdk
import qiime.core.type.grammar as grammar
from qiime.core.type import is_semantic_type

from .data_layout import DataLayout


class Plugin:
    def __init__(self, name, version, website, package):
        self.name = name
        self.version = version
        self.website = website
        self.package = package

        self.methods = PluginMethods(self.name, self.package)
        self.visualizers = PluginVisualizers(self.name)

        self.data_layouts = {}
        self.data_layout_readers = {}
        self.data_layout_writers = {}
        self.types = {}
        self.type_to_data_layouts = {}

    def register_data_layout(self, data_layout):
        if not isinstance(data_layout, DataLayout):
            raise TypeError(
                "`data_layout` must be an instance of type %r, not type %r."
                % (DataLayout.__name__, type(data_layout).__name__))
        # Finalize data layout to prevent additional files from being
        # registered on the data layout after it has been registered.
        data_layout.finalize()
        name = data_layout.name
        version = data_layout.version
        self.data_layouts[(name, version)] = data_layout

    def register_data_layout_reader(self, name, version, view_type, reader):
        self.data_layout_readers[(name, version, view_type)] = reader

    def register_data_layout_writer(self, name, version, view_type, writer):
        self.data_layout_writers[(name, version, view_type)] = writer

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


class PluginMethods(dict):
    def __init__(self, plugin_name, package):
        self._package = package
        self._plugin_name = plugin_name
        super().__init__()

    def register_function(self, function, inputs, parameters, outputs, name,
                          description):
        method = qiime.sdk.Method.from_function(function, inputs, parameters,
                                                outputs, name, description,
                                                plugin_name=self._plugin_name)
        self[method.id] = method

    def register_markdown(self, markdown_filepath):
        markdown_filepath = pkg_resources.resource_filename(
            self._package, markdown_filepath)
        method = qiime.sdk.Method.from_markdown(markdown_filepath,
                                                plugin_name=self._plugin_name)
        self[method.id] = method


class PluginVisualizers(dict):
    def __init__(self, plugin_name):
        self._plugin_name = plugin_name
        super().__init__()

    def register_function(self, function, inputs, parameters, name,
                          description):
        visualizer = qiime.sdk.Visualizer.from_function(
            function, inputs, parameters, name, description,
            plugin_name=self._plugin_name)
        self[visualizer.id] = visualizer
