# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import inspect
import collections
import pkg_resources
import types

import frontmatter

import qiime.sdk
import qiime.core.type.grammar as grammar
from qiime.core.callable import MethodCallable, VisualizerCallable
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

    def register_format(self, format):
        if not issubclass(format, FormatBase):
            raise TypeError("%r is not a valid format." % format)
        if format.__name__ in self.formats.keys():
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
