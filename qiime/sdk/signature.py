# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import qiime.core.type


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
        for input_name, (semantic_type, input_view_type) in \
                inputs.items():
            if not qiime.core.type.is_semantic_type(semantic_type):
                raise TypeError("%r for %r is not a semantic qiime type." %
                                (semantic_type, input_name))

        for param_name, (primitive_type, param_view_type) in \
                parameters.items():
            if not qiime.core.type.is_primitive_type(primitive_type):
                raise TypeError("%r for %r is not a primitive qiime type."
                                % (primitive_type, param_name))

        for output_name, (output_semantic_type, output_view_type) in \
                outputs.items():
            if not qiime.core.type.is_semantic_type(output_semantic_type) and \
               output_semantic_type is not qiime.core.type.Visualization:
                raise TypeError("%r for %r is not a semantic qiime type or "
                                "visualization." %
                                (output_semantic_type, output_name))

        self.inputs = inputs
        self.parameters = parameters
        self.outputs = outputs

    def __call__(self, artifacts, arguments):
        # TODO implement me!
        return self.outputs

    def __eq__(self, other):
        return (
            type(self) is type(other) and
            self.inputs == other.inputs and
            self.parameters == other.parameters and
            self.outputs == other.outputs
        )
