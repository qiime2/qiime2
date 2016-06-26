# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import itertools

import qiime.core.type as qtype

# TODO: Integrate this with workflows
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
            if not qtype.is_semantic_type(semantic_type):
                raise TypeError("%r for %r is not a semantic qiime type." %
                                (semantic_type, input_name))

        for param_name, (primitive_type, _) in parameters.items():
            if not qtype.is_primitive_type(primitive_type):
                raise TypeError("%r for %r is not a primitive qiime type."
                                % (primitive_type, param_name))

        for output_name, (output_semantic_type, _) in outputs.items():
            if not qtype.is_semantic_type(output_semantic_type):
                raise TypeError("%r for %r is not a semantic qiime type."
                                % (output_semantic_type, output_name))
            if not output_semantic_type.is_concrete():
                raise TypeError("%r for %r is not a concrete type."
                                % (output_semantic_type, output_name))

        self.inputs = inputs
        self.parameters = parameters
        self.outputs = outputs

    def decode_parameters(self, **kwargs):
        return {key: exp.decode(arg) for key, arg, exp in
                self._dzip(encoded_arguments, self.parameters)}

    def encode_parameters(self, **kwargs):
        # TODO: Do this.
        pass

    def convert_to_view(self, **kwargs):
        for key, (_, view_type) in self.inputs.items():
            try:
                artifact = kwargs[key]
            except KeyError:
                raise TypeError()

            kwargs[key] = artifact.view(view_type)

        return kwargs

    def type_check(self, **kwargs):
        for key, (type_, _) in itertools.chain(self.inputs.items(),
                                               self.parameters.items()):
            if kwargs[key] not in type_:
                raise TypeError()

    def solve(self, **input_types):
        pass

    def __repr__(self):
        return ", ".join("%s : %r" % entry
                         for entry in itertools.chain(self.inputs.items(),
                                                      self.parameters.items())
                        ) + " -> %r" % (tuple(self.outputs.values()),)

    def __eq__(self, other):
        return (type(self) is type(other) and
                self.name == other.name and
                self.inputs == other.inputs and
                self.parameters == other.parameters and
                self.outputs == other.outputs)
