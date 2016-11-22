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

import qiime.sdk
from .grammar import TypeExpression
from .primitive import is_primitive_type
from .semantic import is_semantic_type
from .visualization import Visualization


# Note: Pipeline doesn't exist yet but it is expected to accept one or more
# input semantic types, zero or more parameters, and produce one or more output
# semantic types or Visualization types.
class PipelineSignature:
    def __init__(self, callable, inputs, parameters, outputs):
        """

        Parameters
        ----------
        callable : callable
            Callable with view type annotations on parameters and return.
        inputs : dict
            Parameter name to semantic type.
        parameters : dict
            Parameter name to primitive type.
        outputs : list of tuple
            Each tuple contains the name of the output (str) and its QIIME
            type.

        """
        inputs, parameters, outputs, defaults = self._parse_annotations(
            callable, inputs, parameters, outputs)

        self._assert_valid_inputs(inputs)
        self._assert_valid_parameters(parameters)
        self._assert_valid_outputs(outputs)
        self._assert_valid_defaults(defaults, parameters)

        self.inputs = inputs
        self.parameters = parameters
        self.outputs = outputs
        self.defaults = defaults

    def _parse_annotations(self, callable, inputs, parameters, outputs):
        sig_params = inspect.signature(callable).parameters
        annotations = {n: p.annotation for n, p in sig_params.items()}

        input_types = collections.OrderedDict()
        param_types = collections.OrderedDict()

        # TODO prevent callable signature from mixing artifacts and primitives.
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
            self._get_return_annotation(callable))

        # TODO make sure these are the same length before zipping
        output_view_types = dict(zip(outputs, output_view_types))

        output_types = collections.OrderedDict()
        for output_name, semantic_type in outputs.items():
            view_type = output_view_types[output_name]
            output_types[output_name] = (semantic_type, view_type)

        return input_types, param_types, output_types, defaults

    def _get_return_annotation(self, callable):
        if 'return' not in callable.__annotations__:
            raise TypeError(
                "Callable %r must have a return type annotation."
                % callable.__name__)

        return callable.__annotations__['return']

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
    def __init__(self, callable, inputs, parameters, outputs):
        super().__init__(callable, inputs, parameters, outputs)

        # Assert all output types are semantic types. The parent class is less
        # strict in its output type requirements.
        for output_name, (output_semantic_type, _) in self.outputs.items():
            if not is_semantic_type(output_semantic_type):
                raise TypeError(
                    "Output %r must be a semantic QIIME type, not %r" %
                    (output_name, output_semantic_type))


class VisualizerSignature(PipelineSignature):
    def __init__(self, callable, inputs, parameters):
        if 'output_dir' in inputs or 'output_dir' in parameters:
            raise TypeError(
                "`output_dir` is a reserved parameter name and cannot be used "
                "in `inputs` or `parameters`")

        callable_parameters = \
            list(inspect.signature(callable).parameters.keys())
        if len(callable_parameters) > 1:
            first_parameter = callable_parameters[0]
            if first_parameter != 'output_dir':
                raise TypeError(
                    "Visualizer callable must have `output_dir` as its first "
                    "argument, not %r" % first_parameter)
        else:
            raise TypeError(
                "Visualizer callable must have at least two arguments, not %d"
                % len(callable_parameters))

        outputs = [('visualization', Visualization)]
        super().__init__(callable, inputs, parameters, outputs)

        _, (_, return_type) = next(iter(self.outputs.items()))
        if return_type is not None:
            raise TypeError(
                "Visualizer callable %r cannot return anything. Its return "
                "annotation must be None, not %r. Write output to "
                "`output_dir`." % (callable.__name__, return_type))
