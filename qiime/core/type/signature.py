# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import itertools

from .primitive import is_primitive_type
from .semantic import is_semantic_type
from .visualization import Visualization


class Signature:
    def __init__(self, inputs, parameters, outputs):
        """
        Parameters
        ----------
        inputs : dict
            Parameter name to tuple of semantic type and view type.
        parameters : dict
            Parameter name to tuple of primitive type and view type.
        outputs : collections.OrderedDict
            Named output to tuple of semantic type and view type.

        """
        for input_name, (semantic_type, _) in inputs.items():
            if not is_semantic_type(semantic_type):
                raise TypeError("%r for %r is not a semantic QIIME type." %
                                (semantic_type, input_name))

        for param_name, (primitive_type, _) in parameters.items():
            if not is_primitive_type(primitive_type):
                raise TypeError("%r for %r is not a primitive QIIME type."
                                % (primitive_type, param_name))

        for output_name, (output_semantic_type, _) in outputs.items():
            if not (is_semantic_type(output_semantic_type) or
                    output_semantic_type == Visualization):
                raise TypeError("%r for %r is not a semantic QIIME type."
                                % (output_semantic_type, output_name))
            if not output_semantic_type.is_concrete():
                raise TypeError("%r for %r is not a concrete type."
                                % (output_semantic_type, output_name))

        self.inputs = inputs
        self.parameters = parameters
        self.outputs = outputs

    def decode_parameters(self, **kwargs):
        params = {}
        for key, (type_expr, _) in self.parameters.items():
            params[key] = type_expr.decode(kwargs[key])
        return params

    def check_types(self, **kwargs):
        for key, (type_, _) in itertools.chain(self.inputs.items(),
                                               self.parameters.items()):
            if kwargs[key] not in type_:
                raise TypeError("Argument to %r is not a subtype of %r."
                                % (key, type_))

    def solve_output(self, **input_types):
        return self.outputs

    def __repr__(self):
        return ", ".join("%s : %r" % entry
                         for entry in itertools.chain(self.inputs.items(),
                                                      self.parameters.items())
                         ) + " -> %r" % (tuple(self.outputs.values()),)

    def __eq__(self, other):
        return (type(self) is type(other) and
                self.inputs == other.inputs and
                self.parameters == other.parameters and
                self.outputs == other.outputs)

    def __ne__(self, other):
        return not (self == other)
