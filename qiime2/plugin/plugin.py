# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import inspect
import types

import qiime2.sdk
import qiime2.core.type.grammar as grammar
from qiime2.core.validate import ValidationObject
from qiime2.plugin.model import DirectoryFormat
from qiime2.plugin.model.base import FormatBase
from qiime2.core.type import is_semantic_type
from qiime2.core.util import get_view_name


TransformerRecord = collections.namedtuple(
    'TransformerRecord', ['transformer', 'plugin', 'citations'])
SemanticTypeRecord = collections.namedtuple(
    'SemanticTypeRecord', ['semantic_type', 'plugin'])
SemanticTypeFragmentRecord = collections.namedtuple(
    'SemanticTypeFragmentRecord', ['fragment', 'plugin'])
FormatRecord = collections.namedtuple('FormatRecord', ['format', 'plugin'])
ViewRecord = collections.namedtuple(
    'ViewRecord', ['name', 'view', 'plugin', 'citations'])
TypeFormatRecord = collections.namedtuple(
    'TypeFormatRecord', ['type_expression', 'format', 'plugin', 'description',
                         'examples'])
ValidatorRecord = collections.namedtuple(
    'ValidatorRecord', ['validator', 'view', 'plugin', 'context'])


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
        self.type_fragments = {}
        self.transformers = {}
        self.type_formats = []
        self.validators = {}

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

    @property
    def types(self):
        types = {}

        for record in self.type_formats:
            for type_ in record.type_expression:
                types[str(type_)] = \
                    SemanticTypeRecord(semantic_type=type_, plugin=self)

        return types

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

    def register_validator(self, semantic_expression):
        if not is_semantic_type(semantic_expression):
            raise TypeError('%s is not a Semantic Type' % semantic_expression)

        def decorator(validator):

            validator_signature = inspect.getfullargspec(validator)

            if 'data' not in validator_signature.annotations:
                raise TypeError('No expected view type provided as annotation'
                                ' for `data` variable in %r.' %
                                (validator.__name__))

            if not ['data', 'level'] == validator_signature.args:
                raise TypeError('The function signature: %r does not contain'
                                ' the required arguments and only the required'
                                ' arguments: %r' % (
                                    validator_signature.args,
                                    ['data', 'level']))

            for semantic_type in semantic_expression:
                if semantic_type not in self.validators:
                    self.validators[semantic_type] = \
                        ValidationObject(semantic_type)

                self.validators[semantic_type].add_validator(
                    ValidatorRecord(
                        validator=validator,
                        view=validator.__annotations__['data'],
                        plugin=self,
                        context=semantic_expression))
            return validator
        return decorator

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

    def register_semantic_types(self, *type_fragments):
        for type_fragment in type_fragments:
            if not is_semantic_type(type_fragment):
                raise TypeError("%r is not a semantic type." % type_fragment)

            if not (isinstance(type_fragment, grammar.IncompleteExp) or
                    (type_fragment.is_concrete() and
                    not type_fragment.fields)):
                raise ValueError("%r is not a semantic type symbol."
                                 % type_fragment)

            if type_fragment.name in self.type_fragments:
                raise ValueError("Duplicate semantic type symbol %r."
                                 % type_fragment)

            self.type_fragments[type_fragment.name] = \
                SemanticTypeFragmentRecord(
                    fragment=type_fragment, plugin=self)

    def _register_artifact_class(self, semantic_type, directory_format,
                                 description, examples):
        if not issubclass(directory_format, DirectoryFormat):
            raise TypeError("%r is not a directory format." % directory_format)
        if not is_semantic_type(semantic_type):
            raise TypeError("%r is not a semantic type." % semantic_type)
        if not is_semantic_type(semantic_type):
            # Evan: this was copied from register_semantic_type_to_format
            # but is unreachable so should be deleted, right?
            raise ValueError("%r is not a semantic type expression."
                             % semantic_type)

        for t in semantic_type:
            if t.predicate is not None:
                raise ValueError("%r has a predicate, differentiating format"
                                 " on predicate is not supported.")

        if description is None:
            description = ""
        if examples is None:
            examples = {}

        self.type_formats.append(TypeFormatRecord(
            type_expression=semantic_type, format=directory_format,
            plugin=self, description=description, examples=examples))

    def register_semantic_type_to_format(self, semantic_type,
                                         artifact_format=None,
                                         directory_format=None):
        # Handle the deprecated parameter name, artifact_format. This is being
        # replaced with directory_format for clarity.
        if artifact_format is not None and directory_format is not None:
            raise ValueError('directory_format and artifact_format were both'
                             'provided when registering artifact class %s.'
                             'Please provide directory_format only as '
                             'artifact_format is deprecated.'
                             % str(semantic_type))
        elif artifact_format is None and directory_format is None:
            raise ValueError('directory_format or artifact_format must be '
                             'provided when registering artifact class %s.'
                             'Please provide directory_format only as '
                             'artifact_format is deprecated.'
                             % str(semantic_type))
        else:
            directory_format = directory_format or artifact_format

        self._register_artifact_class(semantic_type=semantic_type,
                                      directory_format=directory_format,
                                      description=None, examples= None)

    def register_artifact_class(self, semantic_type, directory_format,
                                description=None, examples=None):
        # Evan, better way to determine number of types in the type expression?
        if len(list(semantic_type)) > 1:
            raise TypeError("Only a single type can be registered at a time "
                            "with register_artifact_class. Registration "
                            "attempted for %s." % str(semantic_type))
        self._register_artifact_class(
            semantic_type, directory_format, description, examples)



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
