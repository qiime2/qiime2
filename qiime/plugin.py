# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import inspect
import pkg_resources

import qiime.sdk
import qiime.core.type.grammar as grammar
from qiime.core.type import (SemanticType, is_semantic_type, Int, Str, Float,
                             Visualization)

__all__ = ['load', 'Plugin', 'Int', 'Str', 'Float', 'SemanticType']


def load(plugin_name, plugin_entry_point_name=None):
    # TODO centralize 'qiime.plugin', it is used here and in PluginManager
    plugin_group = 'qiime.plugin'

    if plugin_entry_point_name is None:
        plugin_entry_point_name = plugin_name

    try:
        plugin = pkg_resources.load_entry_point(
            plugin_name, plugin_group, plugin_entry_point_name)
    except ImportError:
        try:
            plugin_entry_map = pkg_resources.get_entry_map(plugin_name)
        except pkg_resources.DistributionNotFoundError:
            raise ImportError("Plugin %r is not installed." % plugin_name)

        if plugin_group not in plugin_entry_map:
            raise ImportError(
                "Plugin %r is not a valid QIIME plugin. A valid QIIME plugin "
                "must define its entry points under the %r entry point group."
                % (plugin_name, plugin_group))
        else:
            plugin_entry_point_names = set(plugin_entry_map[plugin_group])
            raise ImportError(
                "Could not find entry point name %r in plugin %r. Available "
                "entry point names: %r" % (plugin_entry_point_name,
                                           plugin_entry_point_names))
    else:
        if not isinstance(plugin, Plugin):
            raise ImportError(
                "Plugin %r is not a valid QIIME plugin. Expected type %r, "
                "not %r" % (Plugin.__name__, type(plugin).__name__))
        return plugin


class Plugin:
    def __init__(self, name, version, website, package):
        self.package = package
        self.name = name
        self.version = version
        self.website = website
        self.workflows = {}
        self.visualizations = {}
        self.data_layouts = {}
        self.types = {}
        self.type_to_data_layouts = {}
        self.data_layout_readers = {}
        self.data_layout_writers = {}

    def register_workflow(self, workflow):
        fn = pkg_resources.resource_filename(self.package, workflow)
        w = qiime.sdk.Workflow.from_markdown(fn)
        self.workflows[w.id] = w

    def register_function(self, name, function, inputs, parameters, outputs,
                          doc=""):
        # TODO where is the best place to convert outputs as a list of tuples
        # into an OrderedDict?
        outputs = collections.OrderedDict(outputs)
        w = qiime.sdk.Workflow.from_function(function, inputs, parameters,
                                             outputs, name, doc)
        self.workflows[w.id] = w

    def register_visualization(self, name, function, inputs, parameters,
                               doc=""):
        if 'output_dir' in inputs or 'output_dir' in parameters:
            raise TypeError(
                "`output_dir` is a reserved parameter name and cannot be used "
                "in `inputs` or `parameters`")
        parameters['output_dir'] = qiime.core.type.Str

        function_parameters = \
            list(inspect.signature(function).parameters.keys())
        if len(function_parameters) > 0:
            first_parameter = function_parameters[0]
            if first_parameter != 'output_dir':
                raise TypeError(
                    "Visualization function must have `output_dir` as its "
                    "first argument, not %r" % first_parameter)
        else:
            raise TypeError(
                "Visualization function must have at least one argument")

        outputs = collections.OrderedDict([('visualization', Visualization)])
        v = qiime.sdk.Workflow.from_function(function, inputs, parameters,
                                             outputs, name, doc)
        self.visualizations[v.id] = v

    def register_data_layout_reader(self, name, version, view_type, reader):
        self.data_layout_readers[(name, version, view_type)] = reader

    def register_data_layout_writer(self, name, version, view_type, writer):
        self.data_layout_writers[(name, version, view_type)] = writer

    def register_data_layout(self, name, version, validator):
        self.data_layouts[(name, version)] = DataLayout(name, version,
                                                        validator)

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

    def __eq__(self, other):
        return (
            self.package == other.package and
            self.name == other.name and
            self.version == other.version and
            self.website == other.website and
            self.workflows == other.workflows and
            self.data_layouts == other.data_layouts
        )


class DataLayout:
    def __init__(self, name, version, validator):
        # TODO: should this constructor be private since calling it
        # directly would return an unregistered format, and what
        # would you do with that? we're considering this private for
        # now, probably revisit __eq__ if that changes.
        self.name = name
        self.version = version
        self._validator = validator
        self._reader_views = {}
        self._writer_views = {}

    def __eq__(self, other):
        return self.name == other.name and self.version == other.version

    def __ne__(self, other):
        return not self == other

    def validate(self, data_dir):
        return self._validator(data_dir)

    def get_reader_views(self):
        # returns view types that this format can be read into
        return set(self._reader_views)

    def get_writer_views(self):
        # returns view types that this format can be written to
        return set(self._writer_views)

    def reader(self, view_type):
        def decorator(reader_function):
            self._reader_views[view_type] = reader_function
            return reader_function
        return decorator

    @property
    def readers(self):
        return self._reader_views

    def writer(self, view_type):
        def decorator(writer_function):
            self._writer_views[view_type] = writer_function
            return writer_function
        return decorator

    @property
    def writers(self):
        return self._writer_views
