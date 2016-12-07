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
import copy

import qiime.sdk
from .grammar import TypeExpression
from .primitive import is_primitive_type
from .semantic import is_semantic_type
from .visualization import Visualization
from ..util import ImmutableBase


class _NoValue:
    def __repr__(self):
        return "NOVALUE"


class ParameterSpec(ImmutableBase):
    def __init__(self, **kwargs):
        self.NOVALUE = _NoValue()
        self.qiime_type = kwargs.pop('qiime_type', self.NOVALUE)
        self.view_type = kwargs.pop('view_type', self.NOVALUE)
        self.default = kwargs.pop('default', self.NOVALUE)
        if kwargs:
            raise TypeError("Got unexpected keyword argument(s) %r"
                            % list(kwargs))
        self._freeze_()

    def has_qiime_type(self):
        return self.qiime_type is not self.NOVALUE

    def has_view_type(self):
        return self.view_type is not self.NOVALUE

    def has_default(self):
        return self.default is not self.NOVALUE

    def __repr__(self):
        return ("ParameterSpec(qiime_type=%r, view_type=%r, default=%r)"
                % (self.qiime_type, self.view_type, self.default))

    def __eq__(self, other):
        if self.has_qiime_type() != other.has_qiime_type():
            return False
        if self.has_qiime_type() and not self.qiime_type.equals(
                other.qiime_type):
            return False

        if self.has_view_type() != other.has_view_type():
            return False
        if self.has_view_type() and self.view_type is not other.view_type:
            return False

        if self.has_default() != other.has_default():
            return False
        if self.has_default() and self.default != other.default:
            return False

        return True

    def __ne__(self, other):
        return not (self == other)


# Note: Pipeline doesn't exist yet but it is expected to accept one or more
# input semantic types, zero or more parameters, and produce one or more output
# semantic types or Visualization types.
class PipelineSignature:
    provided_args = ('ctx',)

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
        inputs, parameters, outputs = self._parse_signature(
            callable, inputs, parameters, outputs)

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

    def _parse_signature(self, callable, inputs, parameters, outputs):
        # Copy so we can "exhaust" the collections and check for missing params
        inputs = copy.copy(inputs)
        parameters = copy.copy(parameters)

        annotated_inputs = collections.OrderedDict()
        annotated_parameters = collections.OrderedDict()
        annotated_outputs = collections.OrderedDict()

        in_parameter_section = False
        for name, parameter in inspect.signature(callable).parameters.items():
            if (parameter.kind == parameter.VAR_POSITIONAL
                    or parameter.kind == parameter.VAR_KEYWORD):
                raise TypeError("Variadic definitions are unsupported: %r"
                                % name)
            spec_kwargs = {}
            if parameter.annotation is not parameter.empty:
                spec_kwargs['view_type'] = parameter.annotation
            if parameter.default is not parameter.empty:
                spec_kwargs['default'] = parameter.default

            if name in inputs:
                if in_parameter_section:
                    # Mixing "parameters" into the "input" section is not
                    # allowed
                    raise TypeError("Artifact inputs must come before"
                                    " parameters in callable signature.")
                annotated_inputs[name] = ParameterSpec(
                    qiime_type=inputs.pop(name), **spec_kwargs)
            elif name in parameters:
                in_parameter_section = True
                annotated_parameters[name] = ParameterSpec(
                    qiime_type=parameters.pop(name), **spec_kwargs)
            elif name not in self.provided_args:
                raise TypeError("Parameter in callable without QIIME type: %r"
                                % name)
        # we should have popped both of these empty by this point
        if inputs or parameters:
            raise TypeError("Callable does not have parameter(s): %r"
                            % list(inputs) + list(parameters))

        if 'return' in callable.__annotations__:
            output_views = qiime.core.util.tuplize(
                callable.__annotations__['return'])

            if len(output_views) != len(outputs):
                raise TypeError()

            for (name, qiime_type), view_type in zip(outputs, output_views):
                annotated_outputs[name] = ParameterSpec(qiime_type=qiime_type,
                                                        view_type=view_type)
        else:
            for name, qiime_type in outputs:
                annotated_outputs[name] = ParameterSpec(qiime_type=qiime_type)

        return annotated_inputs, annotated_parameters, annotated_outputs


    def _assert_valid_inputs(self, inputs):
        if len(inputs) == 0:
            raise TypeError("%s requires at least one input" %
                            self.__class__.__name__)

        for input_name, spec in inputs.items():
            if not is_semantic_type(spec.qiime_type):
                raise TypeError(
                    "Input %r must be a semantic QIIME type, not %r" %
                    (input_name, spec.qiime_type))

            if not isinstance(spec.qiime_type, TypeExpression):
                raise TypeError(
                    "Input %r must be a complete semantic type expression, "
                    "not %r" % (input_name, spec.qiime_type))

            if spec.has_default():
                raise ValueError("Input %r must not have a default value" %
                                 input_name)

    def _assert_valid_parameters(self, parameters):
        for param_name, spec in parameters.items():
            if not is_primitive_type(spec.qiime_type):
                raise TypeError(
                    "Parameter %r must be a primitive QIIME type, not %r" %
                    (param_name, spec.qiime_type))

            if not isinstance(spec.qiime_type, TypeExpression):
                raise TypeError(
                    "Parameter %r must be a complete primitive type "
                    "expression, not %r" % (param_name, spec.qiime_type))

            if (spec.has_default()
                    and spec.default is not None
                    and spec.default not in spec.qiime_type):
                raise TypeError("Default value for parameter %r is not of "
                                "semantic QIIME type %r or None."
                                % (param_name, spec.qiime_type))

    def _assert_valid_outputs(self, outputs):
        if len(outputs) == 0:
            raise TypeError("%s requires at least one output" %
                            self.__class__.__name__)

        for output_name, spec in outputs.items():
            if not (is_semantic_type(spec.qiime_type) or
                    spec.qiime_type == Visualization):
                raise TypeError(
                    "Output %r must be a semantic QIIME type or "
                    "Visualization, not %r" %
                    (output_name, spec.qiime_type))

            if not isinstance(spec.qiime_type, TypeExpression):
                raise TypeError(
                    "Output %r must be a complete type expression, not %r" %
                    (output_name, spec.qiime_type))

    def decode_parameters(self, **kwargs):
        params = {}
        for key, spec in self.parameters.items():
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
                if not (spec.has_default()
                        and spec.default is kwargs[name] is None):
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
    provided_args = ()

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
    provided_args = ('output_dir',)

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

    def _assert_valid_outputs(self, outputs):
        super()._assert_valid_outputs(outputs)
        output = outputs['visualization']
        if output.has_view_type() and output.view_type is not None:
            raise TypeError(
                "Visualizer callable cannot return anything. Its return "
                "annotation must be None, not %r. Write output to "
                "`output_dir`." % output.view_type)
