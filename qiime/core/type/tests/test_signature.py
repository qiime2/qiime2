# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import collections

from qiime.core.type.signature import PipelineSignature
from qiime.plugin import Int
from qiime.core.testing.type import IntSequence1


class TestPipelineSignature(unittest.TestCase):

    def test_constructor(self):
        signature = PipelineSignature(
            inputs={'input1': (IntSequence1, list),
                    'input2': (IntSequence1, list)},
            parameters={'param1': (Int, int),
                        'param2': (Int, int)},
            outputs=collections.OrderedDict([
                ('output1', (IntSequence1, list)),
                ('output2', (IntSequence1, list))]))

        self.assertEqual(signature.inputs,
                         {'input1': (IntSequence1, list),
                          'input2': (IntSequence1, list)})
        self.assertEqual(signature.parameters,
                         {'param1': (Int, int), 'param2': (Int, int)})
        self.assertEqual(
            signature.outputs,
            collections.OrderedDict([
                ('output1', (IntSequence1, list)),
                ('output2', (IntSequence1, list))]))

    def test_solve_output(self):
        signature = PipelineSignature(
            inputs={'input1': (IntSequence1, list),
                    'input2': (IntSequence1, list)},
            parameters={'param1': (Int, int),
                        'param2': (Int, int)},
            outputs=collections.OrderedDict([
                ('output1', (IntSequence1, list)),
                ('output2', (IntSequence1, list))])
            )
        self.assertEqual(
            signature.solve_output(**{}),
            collections.OrderedDict([
                ('output1', (IntSequence1, list)),
                ('output2', (IntSequence1, list))]))


if __name__ == "__main__":
    unittest.main()
