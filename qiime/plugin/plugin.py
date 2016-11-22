# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import types

import yaml

import qiime.sdk
import qiime.core.type.grammar as grammar
from qiime.plugin.model import DirectoryFormat
from qiime.plugin.model.base import FormatBase
from qiime.core.type import is_semantic_type


TransformerRecord = collections.namedtuple(
    'TransformerRecord', ['transformer', 'restrict', 'plugin'])
SemanticTypeRecord = collections.namedtuple(
    'SemanticTypeRecord', ['semantic_type', 'plugin'])
FormatRecord = collections.namedtuple('FormatRecord', ['format', 'plugin'])
TypeFormatRecord = collections.namedtuple(
    'TypeFormatRecord', ['type_expression', 'format', 'plugin'])


class Plugin:
    @staticmethod
    def yaml_representer(dumper, data):
        items = [
            ('version', data.version),
            ('website', data.website)
        ]
        return dumper.represent_dict(items)

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

        self.formats = {}
        self.types = {}
        self.transformers = {}
        self.type_formats = []

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

    def register_formats(self, *formats):
        for format in formats:
            if not issubclass(format, FormatBase):
                raise TypeError("%r is not a valid format." % format)
            if format.__name__ in self.formats:
                raise NameError("%r is already a registered format." % format)

            self.formats[format.__name__] = FormatRecord(format=format,
                                                         plugin=self)

    def register_transformer(self, _fn=None, *, restrict=None):
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
                transformer=transformer, restrict=restrict, plugin=self)
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

            if not (isinstance(semantic_type, grammar.CompositeType) or
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
        if not isinstance(semantic_type, grammar.TypeExpression):
            raise ValueError("%r is not a semantic type expression."
                             % semantic_type)
        if semantic_type.predicate is not None:
            raise ValueError("%r has a predicate, differentiating format on"
                             " predicate is not supported.")

        self.type_formats.append(TypeFormatRecord(
            type_expression=semantic_type,
            format=artifact_format, plugin=self))


yaml.add_representer(Plugin, Plugin.yaml_representer)


class PluginActions(dict):
    _subpackage = None

    def __init__(self, plugin):
        self._plugin = plugin
        self._package = 'qiime.plugins.%s.%s' % (
            self._plugin.name.replace('-', '_'), self._subpackage)
        super().__init__()


class PluginMethods(PluginActions):
    _subpackage = 'methods'

    # TODO is `register` a better name now that functions are the only accepted
    # source (i.e. markdown support is gone)?
    def register_function(self, function, inputs, parameters, outputs, name,
                          description):
        method = qiime.sdk.Method._init(function, inputs, parameters, outputs,
                                        self._package, name, description)
        self[method.id] = method


class PluginVisualizers(PluginActions):
    _subpackage = 'visualizers'

    def register_function(self, function, inputs, parameters, name,
                          description):
        visualizer = qiime.sdk.Visualizer._init(function, inputs, parameters,
                                                self._package, name,
                                                description)
        self[visualizer.id] = visualizer
