# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
import sys
import importlib.machinery

import qiime.sdk
import qiime.core.data_layout
import qiime.core.type.grammar as grammar
from qiime.core.type import (SemanticType, is_semantic_type, Int, Str, Float,
                             Color, Metadata, MetadataCategory, Properties,
                             Range, Choices, Arguments)

__all__ = ['Plugin', 'Int', 'Str', 'Float', 'Color', 'Metadata',
           'MetadataCategory', 'SemanticType', 'Properties', 'Range',
           'Choices', 'Arguments']
__path__ = []


class Plugin:
    def __init__(self, name, version, website, package):
        self.name = name
        self.version = version
        self.website = website
        self.package = package

        self.methods = PluginMethods(self.name, self.package)
        self.visualizers = PluginVisualizers(self.name)

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
    def __init__(self, plugin, package):
        self._package = package
        self._plugin = plugin
        super().__init__()

    def register_function(self, function, inputs, parameters, outputs, name,
                          description):
        method = qiime.sdk.Method.from_function(function, inputs, parameters,
                                                outputs, name, description,
                                                plugin=self._plugin)
        self[method.id] = method

    def register_markdown(self, markdown_filepath):
        markdown_filepath = pkg_resources.resource_filename(
            self._package, markdown_filepath)
        method = qiime.sdk.Method.from_markdown(markdown_filepath,
                                                plugin=self._plugin)
        self[method.id] = method


class PluginVisualizers(dict):
    def __init__(self, plugin):
        self._plugin = plugin
        super().__init__()

    def register_function(self, function, inputs, parameters, name,
                          description):
        visualizer = qiime.sdk.Visualizer.from_function(
            function, inputs, parameters, name, description,
            plugin=self._plugin)
        self[visualizer.id] = visualizer


class QIIMEArtifactAPIImporter:
    def _plugin_lookup(self, plugin_name):
        import qiime.sdk
        pm = qiime.sdk.PluginManager()
        lookup = {s.replace('-', '_'): s for s in pm.plugins}
        if plugin_name not in lookup:
            return None
        return pm.plugins[lookup[plugin_name]]

    def find_spec(self, name, path=None, target=None):
        # Don't waste time doing anything if it's not a qiime plugin
        if not name.startswith('qiime.plugin.'):
            return None

        if target is not None:
            # TODO: experiment with this to see if it is possible
            raise ImportError("Reloading the QIIME Artifact API is not"
                              " currently supported.")

        # We couldn't care less about path, it is useless to us
        # (It is the __path__ of the parent module)

        fqn = name.split('.')
        plugin_details = fqn[2:]  # fqn[len(['qiime', 'plugin']):]
        plugin_name = plugin_details[0]

        plugin = self._plugin_lookup(plugin_name)
        if plugin is None or len(plugin_details) > 2:
            return None

        if len(plugin_details) == 1:
            return self._make_spec(name, plugin)
        elif plugin_details[1] == 'visualizers':
            return self._make_spec(name, plugin, 'visualizers')
        elif plugin_details[1] == 'methods':
            return self._make_spec(name, plugin, 'methods')
        return None

    def _make_spec(self, name, plugin, action=None):
        # See PEP 451 for explanation of what is happening:
        # https://www.python.org/dev/peps/pep-0451/#modulespec
        return importlib.machinery.ModuleSpec(
            name,
            loader=self,
            origin='generated QIIME API',
            loader_state={'plugin': plugin, 'action': action},
            is_package=action is None
        )

    def exec_module(self, module):
        spec = module.__spec__
        plugin = spec.loader_state['plugin']
        action = spec.loader_state['action']

        if action is None:
            module.methods = importlib.import_module('.methods',
                                                     package=spec.name)
            module.visualizers = importlib.import_module('.visualizers',
                                                         package=spec.name)
        else:
            actions = getattr(plugin, action)
            for key, value in actions.items():
                setattr(module, key, value)

sys.meta_path += [QIIMEArtifactAPIImporter()]
