# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import collections

import qiime.core.testing
from qiime.core.type.signature import Signature
from qiime.plugin import Int


class TestSignature(unittest.TestCase):

    def test_constructor(self):
        signature = Signature(
            inputs={'input1': (qiime.core.testing.TestType, list),
                    'input2': (qiime.core.testing.TestType, list)},
            parameters={'param1': (Int, int),
                        'param2': (Int, int)},
            outputs=collections.OrderedDict([
                ('output1', (qiime.core.testing.TestType, list)),
                ('output2', (qiime.core.testing.TestType, list))]))

        self.assertEqual(signature.inputs,
                         {'input1': (qiime.core.testing.TestType, list),
                          'input2': (qiime.core.testing.TestType, list)})
        self.assertEqual(signature.parameters,
                         {'param1': (Int, int), 'param2': (Int, int)})
        self.assertEqual(
            signature.outputs,
            collections.OrderedDict([
                ('output1', (qiime.core.testing.TestType, list)),
                ('output2', (qiime.core.testing.TestType, list))]))

    def test_solve_output(self):
        signature = Signature(
            inputs={'input1': (qiime.core.testing.TestType, list),
                    'input2': (qiime.core.testing.TestType, list)},
            parameters={'param1': (Int, int),
                        'param2': (Int, int)},
            outputs=collections.OrderedDict([
                ('output1', (qiime.core.testing.TestType, list)),
                ('output2', (qiime.core.testing.TestType, list))])
            )
        self.assertEqual(
            signature.solve_output(**{}),
            collections.OrderedDict([
                ('output1', (qiime.core.testing.TestType, list)),
                ('output2', (qiime.core.testing.TestType, list))]))


if __name__ == "__main__":
    unittest.main()
