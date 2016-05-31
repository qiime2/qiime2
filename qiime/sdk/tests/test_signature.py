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
from qiime.sdk.workflow import Signature
from qiime.plugin import Int


class TestSignature(unittest.TestCase):

    def test_constructor(self):
        signature = Signature(
            'Bogus human-readable name',
            inputs={'input1': (qiime.core.testing.TestType, list),
                    'input2': (qiime.core.testing.TestType, list),
                    'param1': (Int, int),
                    'param2': (Int, int)},
            outputs=collections.OrderedDict([
                ('output1', (qiime.core.testing.TestType, list)),
                ('output2', (qiime.core.testing.TestType, list))]))

        self.assertEqual(signature.name, 'Bogus human-readable name')
        self.assertEqual(signature.input_artifacts,
                         {'input1': (qiime.core.testing.TestType, list),
                          'input2': (qiime.core.testing.TestType, list)})
        self.assertEqual(signature.input_parameters,
                         {'param1': (Int, int), 'param2': (Int, int)})
        self.assertEqual(
            signature.output_artifacts,
            collections.OrderedDict([
                ('output1', (qiime.core.testing.TestType, list)),
                ('output2', (qiime.core.testing.TestType, list))]))

    def test_call(self):
        signature = Signature(
            'Bogus human-readable name',
            inputs={'input1': (qiime.core.testing.TestType, list),
                    'input2': (qiime.core.testing.TestType, list),
                    'param1': (Int, int),
                    'param2': (Int, int)},
            outputs=collections.OrderedDict([
                ('output1', (qiime.core.testing.TestType, list)),
                ('output2', (qiime.core.testing.TestType, list))])
            )
        self.assertEqual(
            signature({}, {}),
            collections.OrderedDict([
                ('output1', (qiime.core.testing.TestType, list)),
                ('output2', (qiime.core.testing.TestType, list))]))


if __name__ == "__main__":
    unittest.main()
