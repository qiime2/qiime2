# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import inspect
import copy
import itertools

import qiime2.sdk
from .grammar import TypeExp, UnionExp
from .meta import TypeVarExp
from .collection import List, Set
from .primitive import infer_primitive_type
from .visualization import Visualization
from . import meta
from .util import is_semantic_type, is_primitive_type
from ..util import ImmutableBase


class __NoValueMeta(type):
    def __repr__(self):
        return "NOVALUE"


# This sentinel is a class so that it retains the correct memory address when
# pickled
class _NOVALUE(metaclass=__NoValueMeta):
    pass


class ParameterSpec(ImmutableBase):
    NOVALUE = _NOVALUE

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

    def duplicate(self, **kwargs):
        qiime_type = kwargs.pop('qiime_type', self.qiime_type)
        view_type = kwargs.pop('view_type', self.view_type)
        default = kwargs.pop('default', self.default)
        description = kwargs.pop('description', self.description)
        if kwargs:
            raise TypeError("Unknown arguments: %r" % kwargs)

        return ParameterSpec(qiime_type, view_type, default, description)

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
        inputs, parameters, outputs, signature_order = \
            self._parse_signature(callable, inputs, parameters, outputs,
                                  input_descriptions, parameter_descriptions,
                                  output_descriptions)

        self._assert_valid_inputs(inputs)
        self._assert_valid_parameters(parameters)
        self._assert_valid_outputs(outputs)
        self._assert_valid_views(inputs, parameters, outputs)

        self.inputs = inputs
        self.parameters = parameters
        self.outputs = outputs
        self.signature_order = signature_order

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
        signature_order = collections.OrderedDict()

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
                description = input_descriptions.pop(name,
                                                     ParameterSpec.NOVALUE)
                param_spec = ParameterSpec(
                    qiime_type=inputs.pop(name), view_type=view_type,
                    default=default, description=description)
                annotated_inputs[name] = param_spec
                signature_order[name] = param_spec
            elif name in parameters:
                description = parameter_descriptions.pop(name,
                                                         ParameterSpec.NOVALUE)
                param_spec = ParameterSpec(
                    qiime_type=parameters.pop(name), view_type=view_type,
                    default=default, description=description)
                annotated_parameters[name] = param_spec
                signature_order[name] = param_spec
            elif name not in self.builtin_args:
                raise TypeError("Parameter in callable without QIIME type:"
                                " %r" % name)
        # we should have popped both of these empty by this point
        if inputs or parameters:
            raise TypeError("Callable does not have parameter(s): %r"
                            % (list(inputs) + list(parameters)))

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

        return (annotated_inputs, annotated_parameters, annotated_outputs,
                signature_order)

    def _assert_valid_inputs(self, inputs):
        for input_name, spec in inputs.items():
            if not is_semantic_type(spec.qiime_type):
                raise TypeError(
                    "Input %r must be a semantic QIIME type, not %r"
                    % (input_name, spec.qiime_type))

            if not isinstance(spec.qiime_type, (TypeExp, UnionExp)):
                raise TypeError(
                    "Input %r must be a complete semantic type expression, "
                    "not %r" % (input_name, spec.qiime_type))

            if spec.has_default() and spec.default is not None:
                raise ValueError(
                    "Input %r has a default value of %r. Only a default "
                    "value of `None` is supported for inputs."
                    % (input_name, spec.default))

            for var_selector in meta.select_variables(spec.qiime_type):
                var = var_selector(spec.qiime_type)
                if not var.input:
                    raise TypeError("An output variable has been associated"
                                    " with an input type: %r"
                                    % spec.qiime_type)

    def _assert_valid_parameters(self, parameters):
        for param_name, spec in parameters.items():
            if not is_primitive_type(spec.qiime_type):
                raise TypeError(
                    "Parameter %r must be a primitive QIIME type, not %r"
                    % (param_name, spec.qiime_type))

            if not isinstance(spec.qiime_type, (TypeExp, UnionExp)):
                raise TypeError(
                    "Parameter %r must be a complete primitive type "
                    "expression, not %r" % (param_name, spec.qiime_type))

            if (spec.has_default() and
                    spec.default is not None and
                    spec.default not in spec.qiime_type):
                raise TypeError("Default value for parameter %r is not of "
                                "semantic QIIME type %r or `None`."
                                % (param_name, spec.qiime_type))

            for var_selector in meta.select_variables(spec.qiime_type):
                var = var_selector(spec.qiime_type)
                if not var.input:
                    raise TypeError("An output variable has been associated"
                                    " with an input type: %r"
                                    % spec.qiime_type)

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

            if not isinstance(spec.qiime_type, (TypeVarExp, TypeExp)):
                raise TypeError(
                    "Output %r must be a complete type expression, not %r"
                    % (output_name, spec.qiime_type))

            for var_selector in meta.select_variables(spec.qiime_type):
                var = var_selector(spec.qiime_type)
                if not var.output:
                    raise TypeError("An input variable has been associated"
                                    " with an input type: %r")

    def _assert_valid_views(self, inputs, parameters, outputs):
        for name, spec in itertools.chain(inputs.items(),
                                          parameters.items(),
                                          outputs.items()):
            if spec.has_view_type():
                raise TypeError(
                    " Pipelines do not support function annotations (found one"
                    " for parameter: %r)." % name)

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
        for name, spec in self.signature_order.items():
            parameter = kwargs[name]
            # A type mismatch is unacceptable unless the value is None
            # and this parameter's default value is None.
            if ((parameter not in spec.qiime_type) and
                    not (spec.has_default() and spec.default is None
                         and parameter is None)):

                if isinstance(parameter, qiime2.sdk.Visualization):
                    raise TypeError(
                        "Parameter %r received a Visualization as an "
                        "argument. Visualizations may not be used as inputs."
                        % name)

                elif isinstance(parameter, qiime2.sdk.Artifact):
                    raise TypeError(
                        "Parameter %r requires an argument of type %r. An "
                        "argument of type %r was passed." % (
                            name, spec.qiime_type, parameter.type))

                elif isinstance(parameter, qiime2.Metadata):
                    raise TypeError(
                        "Parameter %r received Metadata as an "
                        "argument, which is incompatible with parameter "
                        "type: %r" % (name, spec.qiime_type))

                else:  # handle primitive types
                    raise TypeError(
                        "Parameter %r received %r as an argument, which is "
                        "incompatible with parameter type: %r"
                        % (name, parameter, spec.qiime_type))

    def solve_output(self, **kwargs):
        solved_outputs = None
        for _, spec in itertools.chain(self.inputs.items(),
                                       self.parameters.items(),
                                       self.outputs.items()):
            if list(meta.select_variables(spec.qiime_type)):
                break  # a variable exists, do the hard work
        else:
            # no variables
            solved_outputs = self.outputs

        if solved_outputs is None:
            inputs = {**{k: s.qiime_type for k, s in self.inputs.items()},
                      **{k: s.qiime_type for k, s in self.parameters.items()}}
            outputs = {k: s.qiime_type for k, s in self.outputs.items()}
            input_types = {
                k: self._infer_type(k, v) for k, v in kwargs.items()}

            solved = meta.match(input_types, inputs, outputs)
            solved_outputs = collections.OrderedDict(
                (k, s.duplicate(qiime_type=solved[k]))
                for k, s in self.outputs.items())

        for output_name, spec in solved_outputs.items():
            if not spec.qiime_type.is_concrete():
                raise TypeError(
                    "Solved output %r must be a concrete type, not %r" %
                    (output_name, spec.qiime_type))

        return solved_outputs

    def _infer_type(self, key, value):
        if value is None:
            if key in self.inputs:
                return self.inputs[key].qiime_type
            elif key in self.parameters:
                return self.parameters[key].qiime_type
            # Shouldn't happen:
            raise ValueError("Parameter passed not consistent with signature.")
        if type(value) is list:
            inner = UnionExp((self._infer_type(v) for v in value))
            return List[inner.normalize()]
        if type(value) is set:
            inner = UnionExp((self._infer_type(v) for v in value))
            return Set[inner.normalize()]
        if isinstance(value, qiime2.sdk.Artifact):
            return value.type
        else:
            return infer_primitive_type(value)

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
                self.outputs == other.outputs and
                self.signature_order == other.signature_order)

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

    def _assert_valid_views(self, inputs, parameters, outputs):
        for name, spec in itertools.chain(inputs.items(),
                                          parameters.items(),
                                          outputs.items()):
            if not spec.has_view_type():
                raise TypeError("Method is missing a function annotation for"
                                " parameter: %r" % name)


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
                "annotation must be `None`, not %r. Write output to "
                "`output_dir`." % output.view_type)

    def _assert_valid_views(self, inputs, parameters, outputs):
        for name, spec in itertools.chain(inputs.items(), parameters.items()):
            if not spec.has_view_type():
                raise TypeError("Visualizer is missing a function annotation"
                                " for parameter: %r" % name)
