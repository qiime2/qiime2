# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import inspect
import pkg_resources
import types

import frontmatter

import qiime.sdk
import qiime.core.type.grammar as grammar
from qiime.core.callable import MethodCallable, VisualizerCallable
from qiime.core.type import is_semantic_type

from .data_layout import DataLayout


class Plugin:
    def __init__(self, name, version, website, package, citation_text=None,
                 user_support_text=None):
        self.name = name
        self.version = version
        self.website = website
        self.package = package
        if citation_text is None:
            self.citation_text = ('No citation available. Cite plugin '
                                  'website: %s' % self.website)
        else:
            self.citation_text = citation_text
        if user_support_text is None:
            self.user_support_text = ('No user support information available. '
                                      'See plugin website: %s'
                                      % self.website)
        else:
            self.user_support_text = user_support_text

        self.methods = PluginMethods(self)
        self.visualizers = PluginVisualizers(self)

        self.data_layouts = {}
        self.data_layout_readers = {}
        self.data_layout_writers = {}
        self.types = {}
        self.type_to_data_layouts = {}

    @property
    def actions(self):
        # TODO this doesn't handle method/visualizer name collisions. The
        # auto-generated `qiime.plugins.<plugin-name>.actions` API has the same
        # problem. This should be solved at method/visualizer registration
        # time, which will solve the problem for both APIs.
        actions = {}
        actions.update(self.methods)
        actions.update(self.visualizers)
        return types.MappingProxyType(actions)

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


class PluginActions(dict):
    _subpackage = None

    def __init__(self, plugin):
        self._plugin = plugin
        self._package = 'qiime.plugins.%s.%s' % (
            self._plugin.name.replace('-', '_'), self._subpackage)
        super().__init__()

    # Private helpers for use in subclasses:

    def _register_callable(self, callable, name, description, source):
        action = qiime.sdk.Action._from_callable(callable, name, description,
                                                 source)
        self[action.id] = action

    def _get_function_source(self, function):
        try:
            source = inspect.getsource(function)
        except OSError:
            raise TypeError(
                "Cannot retrieve source code for function %r" %
                function.__name__)
        return markdown_source_template % {'source': source}


class PluginMethods(PluginActions):
    _subpackage = 'methods'

    def register_function(self, function, inputs, parameters, outputs, name,
                          description):
        callable = MethodCallable.from_function(function, inputs, parameters,
                                                outputs, self._package)
        source = self._get_function_source(function)

        self._register_callable(callable, name, description, source)

    def register_markdown(self, markdown_filepath):
        markdown_filepath = pkg_resources.resource_filename(
            self._plugin.package, markdown_filepath)

        callable = MethodCallable.from_markdown(markdown_filepath,
                                                self._package)

        with open(markdown_filepath) as fh:
            metadata, source = frontmatter.parse(fh.read())

        self._register_callable(callable, metadata['name'],
                                metadata['description'], source)


class PluginVisualizers(PluginActions):
    _subpackage = 'visualizers'

    def register_function(self, function, inputs, parameters, name,
                          description):
        callable = VisualizerCallable.from_function(function, inputs,
                                                    parameters, self._package)
        source = self._get_function_source(function)

        self._register_callable(callable, name, description, source)


markdown_source_template = """
```python
%(source)s
```
"""
