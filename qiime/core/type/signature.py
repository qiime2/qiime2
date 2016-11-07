# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import inspect
import itertools
import pydoc

import frontmatter

import qiime.sdk
from .grammar import TypeExpression
from .primitive import is_primitive_type
from .semantic import is_semantic_type
from .visualization import Visualization


# Note: Pipeline doesn't exist yet but it is expected to accept one or more
# input semantic types, zero or more parameters, and produce one or more output
# semantic types or Visualization types.
class PipelineSignature:
    def __init__(self, inputs, parameters, defaults, outputs):
        """
        Parameters
        ----------
        inputs : dict
            Parameter name to tuple of semantic type and view type.
        parameters : dict
            Parameter name to tuple of primitive type and view type.
        defaults : dict
            Parameter name to default argument.
        outputs : collections.OrderedDict
            Named output to tuple of QIIME type and view type.
        """
        # TODO ensure `inputs`, `parameters`, and `outputs` are OrderedDicts.
        self._assert_valid_inputs(inputs)
        self._assert_valid_parameters(parameters)
        self._assert_valid_defaults(defaults, parameters)
        self._assert_valid_outputs(outputs)

        self.inputs = inputs
        self.parameters = parameters
        self.outputs = outputs
        self.defaults = defaults

    def _assert_valid_inputs(self, inputs):
        if len(inputs) == 0:
            raise TypeError("%s requires at least one input" %
                            self.__class__.__name__)

        for input_name, (semantic_type, _) in inputs.items():
            if not is_semantic_type(semantic_type):
                raise TypeError(
                    "Input %r must be a semantic QIIME type, not %r" %
                    (input_name, semantic_type))

            if not isinstance(semantic_type, TypeExpression):
                raise TypeError(
                    "Input %r must be a complete semantic type expression, "
                    "not %r" % (input_name, semantic_type))

    def _assert_valid_parameters(self, parameters):
        for param_name, (primitive_type, _) in parameters.items():
            if not is_primitive_type(primitive_type):
                raise TypeError(
                    "Parameter %r must be a primitive QIIME type, not %r" %
                    (param_name, primitive_type))

            if not isinstance(primitive_type, TypeExpression):
                raise TypeError(
                    "Parameter %r must be a complete primitive type "
                    "expression, not %r" % (param_name, primitive_type))

    def _assert_valid_defaults(self, defaults, parameters):
        # only parameters can have defaults and the type must match or be None
        for param_name, default in defaults.items():
            if param_name not in parameters:
                raise ValueError("Input %r must not have a default value" %
                                 param_name)

            if default is not None and \
               default not in parameters[param_name][0]:
                raise TypeError("Default value for parameter %r is not of "
                                "semantic QIIME type %r or None." %
                                (param_name, parameters[param_name][0]))

    def _assert_valid_outputs(self, outputs):
        if len(outputs) == 0:
            raise TypeError("%s requires at least one output" %
                            self.__class__.__name__)

        for output_name, (output_semantic_type, _) in outputs.items():
            if not (is_semantic_type(output_semantic_type) or
                    output_semantic_type == Visualization):
                raise TypeError(
                    "Output %r must be a semantic QIIME type or "
                    "Visualization, not %r" %
                    (output_name, output_semantic_type))

            if not isinstance(output_semantic_type, TypeExpression):
                raise TypeError(
                    "Output %r must be a complete type expression, not %r" %
                    (output_name, output_semantic_type))

    def decode_parameters(self, **kwargs):
        params = {}
        for key, (type_expr, _) in self.parameters.items():
            params[key] = type_expr.decode(kwargs[key])
        return params

    def check_types(self, **kwargs):
        for key, (type_, _) in self.inputs.items():
            if kwargs[key] not in type_:
                raise TypeError("Argument to input %r is not a subtype of"
                                " %r." % (key, type_))

        for key, (type_, _) in self.parameters.items():
            if kwargs[key] not in type_:
                # A type mismatch is unacceptable unless the value is None
                # and this parameter's default value is None.
                if not (key in self.defaults and
                        self.defaults[key] is None and
                        kwargs[key] is None):
                    raise TypeError("Argument to parameter %r is not a "
                                    "subtype of %r." % (key, type_))

    def solve_output(self, **input_types):
        # TODO implement solving here. The check for concrete output types may
        # be unnecessary here if the signature's constructor can detect
        # unsolvable signatures and ensure that solving will always produce
        # concrete output types.
        solved_outputs = self.outputs

        for output_name, (output_semantic_type, _) in solved_outputs.items():
            if not output_semantic_type.is_concrete():
                raise TypeError(
                    "Solved output %r must be a concrete type, not %r" %
                    (output_name, output_semantic_type))

        return solved_outputs

    def __repr__(self):
        inputs = []
        for name, type in itertools.chain(self.inputs.items(),
                                          self.parameters.items()):
            inputs.append("%s : %r" % (name, type))
        return ", ".join(inputs) + " -> %r" % (tuple(self.outputs.values()),)

    def __eq__(self, other):
        return (type(self) is type(other) and
                self.inputs == other.inputs and
                self.parameters == other.parameters and
                self.defaults == other.defaults and
                self.outputs == other.outputs)

    def __ne__(self, other):
        return not (self == other)


class MethodSignature(PipelineSignature):
    @classmethod
    def from_function(cls, function, inputs, parameters, outputs):
        return function_to_signature(
            cls, function, inputs, parameters, outputs)

    @classmethod
    def from_markdown(cls, markdown_filepath):
        with open(markdown_filepath) as fh:
            metadata, _ = frontmatter.parse(fh.read())

        input_types = collections.OrderedDict()
        for input_ in metadata['inputs']:
            name, type_tuple = list(input_.items())[0]
            if len(type_tuple) != 2:
                raise TypeError(
                    "Bad input section formatting, need two items per key"
                )
            input_types[name] = cls._split_type_tuple(type_tuple, 'semantic')

        param_types = collections.OrderedDict()
        for parameter in metadata['parameters']:
            name, type_tuple = list(parameter.items())[0]
            if len(type_tuple) != 2:
                raise TypeError(
                    "Bad parameters section formatting, need two items per key"
                )
            param_types[name] = cls._split_type_tuple(type_tuple, 'primitive')

        output_types = collections.OrderedDict()
        for output in metadata['outputs']:
            name, type_tuple = list(output.items())[0]
            if len(type_tuple) != 2:
                raise TypeError(
                    "Bad outputs section formatting, need two items per key"
                )
            output_types[name] = cls._split_type_tuple(type_tuple, 'semantic')

        # TODO come up with a nice way to format default values in markdown
        return cls(input_types, param_types, {}, output_types)

    @classmethod
    def _split_type_tuple(cls, type_tuple, expect):
        # Split the type expression into its components: the semantic_type
        # and the view_type.
        semantic_type, view_type = type_tuple
        view_type = cls._parse_view_type(view_type)
        semantic_type = qiime.sdk.parse_type(semantic_type, expect=expect)
        return semantic_type, view_type

    # TODO this isn't safe nor robust, and not how we want to ultimately handle
    # importing a view type. Refactor to use a similar import mechanism as
    # semantic types when that part is finalized. Hack taken from:
    #     http://stackoverflow.com/a/29831586/3776794
    @classmethod
    def _parse_view_type(cls, view_type_str):
        view_type = pydoc.locate(view_type_str)
        if view_type is None:
            raise ImportError("Could not import view type %r" % view_type_str)
        return view_type

    def __init__(self, inputs, parameters, defaults, outputs):
        super().__init__(inputs, parameters, defaults, outputs)

        # Assert all output types are semantic types. The parent class is less
        # strict in its output type requirements.
        for output_name, (output_semantic_type, _) in self.outputs.items():
            if not is_semantic_type(output_semantic_type):
                raise TypeError(
                    "Output %r must be a semantic QIIME type, not %r" %
                    (output_name, output_semantic_type))


class VisualizerSignature(PipelineSignature):
    @classmethod
    def from_function(cls, function, inputs, parameters):
        if 'output_dir' in inputs or 'output_dir' in parameters:
            raise TypeError(
                "`output_dir` is a reserved parameter name and cannot be used "
                "in `inputs` or `parameters`")

        function_parameters = \
            list(inspect.signature(function).parameters.keys())
        if len(function_parameters) > 1:
            first_parameter = function_parameters[0]
            if first_parameter != 'output_dir':
                raise TypeError(
                    "Visualizer function must have `output_dir` as its first "
                    "argument, not %r" % first_parameter)
        else:
            raise TypeError(
                "Visualizer function must have at least two arguments, not %d"
                % len(function_parameters))

        return_annotation = get_return_annotation(function)
        if return_annotation is not None:
            raise TypeError(
                "Visualizer function %r cannot return anything. Its return "
                "annotation must be None, not %r. Write output to "
                "`output_dir`." % (function.__name__, return_annotation))

        outputs = [('visualization', Visualization)]
        return function_to_signature(cls, function, inputs, parameters,
                                     outputs)

    def __init__(self, inputs, parameters, defaults, outputs):
        super().__init__(inputs, parameters, defaults, outputs)

        # Assert there is exactly one output that is a Visualization type. The
        # parent class is less strict in its output requirements.
        if len(self.outputs) != 1:
            raise TypeError("%s requires exactly one output, not %d" %
                            (self.__class__.__name__, len(self.outputs)))

        output_name, (output_semantic_type, _) = \
            next(iter(self.outputs.items()))
        if output_semantic_type != Visualization:
            raise TypeError(
                "Output %r must be a Visualization QIIME type, not %r" %
                (output_name, output_semantic_type))


# Utility function for parsing a function into a signature. Used to reduce
# duplicate code between MethodSignature.from_function and
# VisualizerSignature.from_function.
def function_to_signature(cls, function, inputs, parameters, outputs):
    sig_params = inspect.signature(function).parameters
    annotations = {n: p.annotation for n, p in sig_params.items()}

    input_types = collections.OrderedDict()
    param_types = collections.OrderedDict()

    # TODO prevent function signature from mixing artifacts and primitives.
    # Artifacts come first followed by primitives.
    for name in sig_params:
        if name in inputs:
            input_types[name] = (inputs[name], annotations[name])
        if name in parameters:
            param_types[name] = (parameters[name], annotations[name])

    defaults = {n: p.default for n, p in sig_params.items()
                if p.default is not p.empty and n in parameters}

    outputs = collections.OrderedDict(outputs)
    output_view_types = qiime.core.util.tuplize(
        get_return_annotation(function))

    # TODO make sure these are the same length before zipping
    output_view_types = dict(zip(outputs, output_view_types))

    output_types = collections.OrderedDict()
    for output_name, semantic_type in outputs.items():
        view_type = output_view_types.get(output_name)
        
        if output_name != "visualization" and view_type is None:
            raise TypeError("Function %r does not describe %r in its return "
                            "annotation." % (function.__name__, output_name))

        output_types[output_name] = (semantic_type, view_type)

    return cls(input_types, param_types, defaults, output_types)


def get_return_annotation(function):
    if 'return' not in function.__annotations__:
        raise TypeError(
            "Function %r must have a return type annotation."
            % function.__name__)

    return function.__annotations__['return']
