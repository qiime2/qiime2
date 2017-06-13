# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import inspect
import copy

import qiime2.sdk
from .grammar import TypeExpression
from .primitive import is_primitive_type
from .semantic import is_semantic_type
from .visualization import Visualization
from ..util import ImmutableBase


class _NoValue:
    def __repr__(self):
        return "NOVALUE"


class ParameterSpec(ImmutableBase):
    NOVALUE = _NoValue()

    def __init__(self, qiime_type=NOVALUE, view_type=NOVALUE, default=NOVALUE,
                 description=NOVALUE):
        self.qiime_type = qiime_type
        self.view_type = view_type
        self.default = default
        self.description = description

        self._freeze_()

    def has_qiime_type(self):
        return self.qiime_type is not self.NOVALUE

    def has_view_type(self):
        return self.view_type is not self.NOVALUE

    def has_default(self):
        return self.default is not self.NOVALUE

    def has_description(self):
        return self.description is not self.NOVALUE

    def __repr__(self):
        return ("ParameterSpec(qiime_type=%r, view_type=%r, default=%r, "
                "description=%r)" % (self.qiime_type, self.view_type,
                                     self.default, self.description))

    def __eq__(self, other):
        return (self.qiime_type == other.qiime_type and
                self.view_type == other.view_type and
                self.default == other.default and
                self.description == other.description)

    def __ne__(self, other):
        return not (self == other)


# Note: Pipeline doesn't exist yet but it is expected to accept zero or more
# input semantic types, zero or more parameters, and produce one or more output
# semantic types or Visualization types.
class PipelineSignature:
    builtin_args = ('ctx',)

    def __init__(self, callable, inputs, parameters, outputs,
                 input_descriptions=None, parameter_descriptions=None,
                 output_descriptions=None):
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
        input_descriptions : dict, optional
            Input name to description string.
        parameter_descriptions : dict, optional
            Parameter name to description string.
        output_descriptions : dict, optional
            Output name to description string.

        """
        inputs, parameters, outputs = self._parse_signature(
            callable, inputs, parameters, outputs, input_descriptions,
            parameter_descriptions, output_descriptions)

        self._assert_valid_inputs(inputs)
        self._assert_valid_parameters(parameters)
        self._assert_valid_outputs(outputs)

        self.inputs = inputs
        self.parameters = parameters
        self.outputs = outputs

    @property
    def defaults(self):
        return collections.OrderedDict([
            (name, spec.default) for name, spec in self.parameters.items()
            if spec.has_default()])

    def _parse_signature(self, callable, inputs, parameters, outputs,
                         input_descriptions=None, parameter_descriptions=None,
                         output_descriptions=None):
        #  Initialize dictionaries if non-existant.
        if input_descriptions is None:
            input_descriptions = {}
        if parameter_descriptions is None:
            parameter_descriptions = {}
        if output_descriptions is None:
            output_descriptions = {}

        # Copy so we can "exhaust" the collections and check for missing params
        inputs = copy.copy(inputs)
        parameters = copy.copy(parameters)
        input_descriptions = copy.copy(input_descriptions)
        parameter_descriptions = copy.copy(parameter_descriptions)
        output_descriptions = copy.copy(output_descriptions)
        builtin_args = list(self.builtin_args)

        annotated_inputs = collections.OrderedDict()
        annotated_parameters = collections.OrderedDict()
        annotated_outputs = collections.OrderedDict()

        in_parameter_section = False
        for name, parameter in inspect.signature(callable).parameters.items():
            if (parameter.kind == parameter.VAR_POSITIONAL or
                    parameter.kind == parameter.VAR_KEYWORD):
                raise TypeError("Variadic definitions are unsupported: %r" %
                                name)

            if builtin_args:
                if builtin_args[0] != name:
                    raise TypeError("Missing builtin argument %r, got %r" %
                                    (builtin_args[0], name))
                builtin_args = builtin_args[1:]
                continue

            view_type = ParameterSpec.NOVALUE
            if parameter.annotation is not parameter.empty:
                view_type = parameter.annotation
            default = ParameterSpec.NOVALUE
            if parameter.default is not parameter.empty:
                default = parameter.default

            if name in inputs:
                if in_parameter_section:
                    # Mixing "parameters" into the "input" section is not
                    # allowed
                    raise TypeError("Artifact inputs must come before"
                                    " parameters in callable signature.")
                description = input_descriptions.pop(name,
                                                     ParameterSpec.NOVALUE)
                annotated_inputs[name] = ParameterSpec(
                    qiime_type=inputs.pop(name), view_type=view_type,
                    default=default, description=description)
            elif name in parameters:
                in_parameter_section = True
                description = parameter_descriptions.pop(name,
                                                         ParameterSpec.NOVALUE)
                annotated_parameters[name] = ParameterSpec(
                    qiime_type=parameters.pop(name), view_type=view_type,
                    default=default, description=description)
            elif name not in self.builtin_args:
                raise TypeError("Parameter in callable without QIIME type:"
                                " %r" % name)
        # we should have popped both of these empty by this point
        if inputs or parameters:
            raise TypeError("Callable does not have parameter(s): %r"
                            % list(inputs) + list(parameters))

        if 'return' in callable.__annotations__:
            output_views = qiime2.core.util.tuplize(
                callable.__annotations__['return'])

            if len(output_views) != len(outputs):
                raise TypeError("Number of registered outputs (%r) does not"
                                " match annotation (%r)" %
                                (len(outputs), len(output_views)))

            for (name, qiime_type), view_type in zip(outputs, output_views):
                description = output_descriptions.pop(name,
                                                      ParameterSpec.NOVALUE)
                annotated_outputs[name] = ParameterSpec(
                    qiime_type=qiime_type, view_type=view_type,
                    description=description)
        else:
            for name, qiime_type in outputs:
                description = output_descriptions.pop(name,
                                                      ParameterSpec.NOVALUE)
                annotated_outputs[name] = ParameterSpec(
                    qiime_type=qiime_type, description=description)

        # we should have popped the descriptions empty by this point
        if input_descriptions or parameter_descriptions or output_descriptions:
            raise TypeError(
                "Callable does not have parameter(s)/output(s) found in "
                "descriptions: %r" % [*input_descriptions,
                                      *parameter_descriptions,
                                      *output_descriptions])

        return annotated_inputs, annotated_parameters, annotated_outputs

    def _assert_valid_inputs(self, inputs):
        for input_name, spec in inputs.items():
            if not is_semantic_type(spec.qiime_type):
                raise TypeError(
                    "Input %r must be a semantic QIIME type, not %r"
                    % (input_name, spec.qiime_type))

            if not isinstance(spec.qiime_type, TypeExpression):
                raise TypeError(
                    "Input %r must be a complete semantic type expression, "
                    "not %r" % (input_name, spec.qiime_type))

            if spec.has_default():
                raise ValueError("Input %r must not have a default value"
                                 % input_name)

    def _assert_valid_parameters(self, parameters):
        for param_name, spec in parameters.items():
            if not is_primitive_type(spec.qiime_type):
                raise TypeError(
                    "Parameter %r must be a primitive QIIME type, not %r"
                    % (param_name, spec.qiime_type))

            if not isinstance(spec.qiime_type, TypeExpression):
                raise TypeError(
                    "Parameter %r must be a complete primitive type "
                    "expression, not %r" % (param_name, spec.qiime_type))

            if (spec.has_default() and
                    spec.default is not None and
                    spec.default not in spec.qiime_type):
                raise TypeError("Default value for parameter %r is not of "
                                "semantic QIIME type %r or None."
                                % (param_name, spec.qiime_type))

    def _assert_valid_outputs(self, outputs):
        if len(outputs) == 0:
            raise TypeError("%s requires at least one output"
                            % self.__class__.__name__)

        for output_name, spec in outputs.items():
            if not (is_semantic_type(spec.qiime_type) or
                    spec.qiime_type == Visualization):
                raise TypeError(
                    "Output %r must be a semantic QIIME type or "
                    "Visualization, not %r"
                    % (output_name, spec.qiime_type))

            if not isinstance(spec.qiime_type, TypeExpression):
                raise TypeError(
                    "Output %r must be a complete type expression, not %r"
                    % (output_name, spec.qiime_type))

    def decode_parameters(self, **kwargs):
        params = {}
        for key, spec in self.parameters.items():
            if (spec.has_default() and
                    spec.default is None and
                    kwargs[key] is None):
                params[key] = None
            else:
                params[key] = spec.qiime_type.decode(kwargs[key])
        return params

    def check_types(self, **kwargs):
        for name, spec in self.inputs.items():
            if kwargs[name] not in spec.qiime_type:
                raise TypeError("Argument to input %r is not a subtype of"
                                " %r." % (name, spec.qiime_type))

        for name, spec in self.parameters.items():
            if kwargs[name] not in spec.qiime_type:
                # A type mismatch is unacceptable unless the value is None
                # and this parameter's default value is None.
                if not (spec.has_default() and
                        spec.default is None and
                        kwargs[name] is None):
                    raise TypeError("Argument to parameter %r is not a "
                                    "subtype of %r." % (name, spec.qiime_type))

    def solve_output(self, **input_types):
        # TODO implement solving here. The check for concrete output types may
        # be unnecessary here if the signature's constructor can detect
        # unsolvable signatures and ensure that solving will always produce
        # concrete output types.
        solved_outputs = self.outputs

        for output_name, spec in solved_outputs.items():
            if not spec.qiime_type.is_concrete():
                raise TypeError(
                    "Solved output %r must be a concrete type, not %r" %
                    (output_name, spec.qiime_type))

        return solved_outputs

    def __repr__(self):
        lines = []
        for group in 'inputs', 'parameters', 'outputs':
            lookup = getattr(self, group)
            lines.append('%s:' % group)
            for name, spec in lookup.items():
                lines.append('    %s: %r' % (name, spec))
        return '\n'.join(lines)

    def __eq__(self, other):
        return (type(self) is type(other) and
                self.inputs == other.inputs and
                self.parameters == other.parameters and
                self.outputs == other.outputs)

    def __ne__(self, other):
        return not (self == other)


class MethodSignature(PipelineSignature):
    builtin_args = ()

    def _assert_valid_outputs(self, outputs):
        super()._assert_valid_outputs(outputs)
        # Assert all output types are semantic types. The parent class is less
        # strict in its output type requirements.
        for output_name, spec in outputs.items():
            if not is_semantic_type(spec.qiime_type):
                raise TypeError(
                    "Output %r must be a semantic QIIME type, not %r" %
                    (output_name, spec.qiime_type))


class VisualizerSignature(PipelineSignature):
    builtin_args = ('output_dir',)

    def __init__(self, callable, inputs, parameters, input_descriptions=None,
                 parameter_descriptions=None):
        outputs = [('visualization', Visualization)]
        output_descriptions = None
        super().__init__(callable, inputs, parameters, outputs,
                         input_descriptions, parameter_descriptions,
                         output_descriptions)

    def _assert_valid_outputs(self, outputs):
        super()._assert_valid_outputs(outputs)
        output = outputs['visualization']
        if output.has_view_type() and output.view_type is not None:
            raise TypeError(
                "Visualizer callable cannot return anything. Its return "
                "annotation must be None, not %r. Write output to "
                "`output_dir`." % output.view_type)
