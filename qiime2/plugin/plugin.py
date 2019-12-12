# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import types

import qiime2.sdk
import qiime2.core.type.grammar as grammar
from qiime2.plugin.model import DirectoryFormat
from qiime2.plugin.model.base import FormatBase
from qiime2.core.type import is_semantic_type
from qiime2.core.util import get_view_name


TransformerRecord = collections.namedtuple(
    'TransformerRecord', ['transformer', 'plugin', 'citations'])
SemanticTypeRecord = collections.namedtuple(
    'SemanticTypeRecord', ['semantic_type', 'plugin'])
FormatRecord = collections.namedtuple('FormatRecord', ['format', 'plugin'])
ViewRecord = collections.namedtuple(
    'ViewRecord', ['name', 'view', 'plugin', 'citations'])
TypeFormatRecord = collections.namedtuple(
    'TypeFormatRecord', ['type_expression', 'format', 'plugin'])


class Plugin:
    def __init__(self, name, version, website, package=None, project_name=None,
                 citation_text=None, user_support_text=None,
                 short_description=None, description=None, citations=None):
        self.id = name.replace('-', '_')
        self.name = name
        self.version = version
        self.website = website

        # Filled in by the PluginManager if not provided.
        self.package = package
        self.project_name = project_name

        if user_support_text is None:
            self.user_support_text = ('Please post to the QIIME 2 forum for '
                                      'help with this plugin: https://forum.'
                                      'qiime2.org')
        else:
            self.user_support_text = user_support_text

        if short_description is None:
            self.short_description = ''
        else:
            self.short_description = short_description

        if description is None:
            self.description = ('No description available. '
                                'See plugin website: %s'
                                % self.website)
        else:
            self.description = description

        if citations is None:
            self.citations = ()
        else:
            self.citations = tuple(citations)

        self.methods = PluginMethods(self)
        self.visualizers = PluginVisualizers(self)
        self.pipelines = PluginPipelines(self)

        self.formats = {}
        self.views = {}
        self.types = {}
        self.transformers = {}
        self.type_formats = []

    def freeze(self):
        pass

    @property
    def actions(self):
        # TODO this doesn't handle method/visualizer name collisions. The
        # auto-generated `qiime2.plugins.<plugin-name>.actions` API has the
        # same problem. This should be solved at method/visualizer registration
        # time, which will solve the problem for both APIs.
        actions = {}
        actions.update(self.methods)
        actions.update(self.visualizers)
        actions.update(self.pipelines)
        return types.MappingProxyType(actions)

    def register_formats(self, *formats, citations=None):
        for format in formats:
            if not issubclass(format, FormatBase):
                raise TypeError("%r is not a valid format." % format)

        self.register_views(*formats, citations=citations)

    def register_views(self, *views, citations=None):
        if citations is None:
            citations = ()
        else:
            citations = tuple(citations)

        for view in views:
            if not isinstance(view, type):
                raise TypeError("%r should be a class." % view)

            is_format = False
            if issubclass(view, FormatBase):
                is_format = True

            name = get_view_name(view)
            if name in self.views:
                raise NameError("View %r is already registered by this "
                                "plugin." % name)

            self.views[name] = ViewRecord(
                name=name, view=view, plugin=self, citations=citations)

            if is_format:
                self.formats[name] = FormatRecord(format=view, plugin=self)

    def register_transformer(self, _fn=None, *, citations=None):
        """
        A transformer has the type Callable[[type], type]
        """
        # `_fn` allows us to figure out if we are called with or without
        # arguments in order to support both:
        # ```
        # @plugin.register_transformer
        # def _(x: A) -> B:
        #    ...
        # ```
        # and
        # ```
        # @plugin.register_transformer(restrict=True)
        # def _(x: A) -> B:
        #   ...
        # ```
        if citations is None:
            citations = ()
        else:
            citations = tuple(citations)

        def decorator(transformer):
            annotations = transformer.__annotations__.copy()
            if len(annotations) != 2:
                raise TypeError("A transformer must only have a single input"
                                " and output annotation.")
            try:
                output = annotations.pop('return')
            except KeyError:
                raise TypeError("A transformer must provide a return type.")

            if type(output) is tuple:
                raise TypeError("A transformer can only return a single type,"
                                " not %r." % (output,))

            input = list(annotations.values())[0]
            if (input, output) in self.transformers:
                raise TypeError("Duplicate transformer (%r) from %r to %r."
                                % (transformer, input, output))
            if input == output:
                raise TypeError("Plugins should not register identity"
                                " transformations (%r, %r to %r)."
                                % (transformer, input, output))

            self.transformers[input, output] = TransformerRecord(
                transformer=transformer, plugin=self, citations=citations)
            return transformer

        if _fn is None:
            return decorator
        else:
            # Apply the decorator as we were applied with a single function
            return decorator(_fn)

    def register_semantic_types(self, *semantic_types):
        for semantic_type in semantic_types:
            if not is_semantic_type(semantic_type):
                raise TypeError("%r is not a semantic type." % semantic_type)

            if not (isinstance(semantic_type, grammar.IncompleteExp) or
                    (semantic_type.is_concrete() and
                    not semantic_type.fields)):
                raise ValueError("%r is not a semantic type symbol."
                                 % semantic_type)

            if semantic_type.name in self.types:
                raise ValueError("Duplicate semantic type symbol %r."
                                 % semantic_type)

            self.types[semantic_type.name] = SemanticTypeRecord(
                semantic_type=semantic_type, plugin=self)

    def register_semantic_type_to_format(self, semantic_type, artifact_format):
        if not issubclass(artifact_format, DirectoryFormat):
            raise TypeError("%r is not a directory format." % artifact_format)
        if not is_semantic_type(semantic_type):
            raise TypeError("%r is not a semantic type." % semantic_type)
        if not is_semantic_type(semantic_type):
            raise ValueError("%r is not a semantic type expression."
                             % semantic_type)
        for t in semantic_type:
            if t.predicate is not None:
                raise ValueError("%r has a predicate, differentiating format"
                                 " on predicate is not supported.")

        self.type_formats.append(TypeFormatRecord(
            type_expression=semantic_type,
            format=artifact_format, plugin=self))


class PluginActions(dict):
    def __init__(self, plugin):
        self._plugin_id = plugin.id
        super().__init__()


class PluginMethods(PluginActions):
    def register_function(self, function, inputs, parameters, outputs, name,
                          description, input_descriptions=None,
                          parameter_descriptions=None,
                          output_descriptions=None, citations=None,
                          deprecated=False, examples=None):
        if citations is None:
            citations = ()
        else:
            citations = tuple(citations)

        if examples is None:
            examples = {}

        method = qiime2.sdk.Method._init(function, inputs, parameters, outputs,
                                         self._plugin_id, name, description,
                                         input_descriptions,
                                         parameter_descriptions,
                                         output_descriptions, citations,
                                         deprecated, examples)
        self[method.id] = method


class PluginVisualizers(PluginActions):
    def register_function(self, function, inputs, parameters, name,
                          description, input_descriptions=None,
                          parameter_descriptions=None, citations=None,
                          deprecated=False, examples=None):
        if citations is None:
            citations = ()
        else:
            citations = tuple(citations)

        if examples is None:
            examples = {}

        visualizer = qiime2.sdk.Visualizer._init(function, inputs, parameters,
                                                 self._plugin_id, name,
                                                 description,
                                                 input_descriptions,
                                                 parameter_descriptions,
                                                 citations, deprecated,
                                                 examples)
        self[visualizer.id] = visualizer


class PluginPipelines(PluginActions):
    def register_function(self, function, inputs, parameters, outputs, name,
                          description, input_descriptions=None,
                          parameter_descriptions=None,
                          output_descriptions=None, citations=None,
                          deprecated=False, examples=None):
        if citations is None:
            citations = ()
        else:
            citations = tuple(citations)

        if examples is None:
            examples = {}

        pipeline = qiime2.sdk.Pipeline._init(function, inputs, parameters,
                                             outputs, self._plugin_id, name,
                                             description, input_descriptions,
                                             parameter_descriptions,
                                             output_descriptions, citations,
                                             deprecated, examples)
        self[pipeline.id] = pipeline
