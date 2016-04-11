# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import collections

from qiime.sdk.workflow import Signature
from qiime.plugin import Int, Type


class DummyType(Type, variant_of=Type.Artifact):
    def load(self, data_reader):
        fh = data_reader.get_file('data.txt')
        model = []
        for line in fh:
            model.append(int(line.rstrip()))
        return model

    def save(self, data, data_writer):
        fh = data_writer.create_file('data.txt')
        for num in data:
            fh.write('%d\n' % num)


class TestSignature(unittest.TestCase):

    def test_constructor(self):
        signature = Signature(
            'Bogus human-readable name',
            inputs={'input1': DummyType,
                    'input2': DummyType,
                    'param1': Int,
                    'param2': Int},
            outputs=collections.OrderedDict([('output1', DummyType),
                                             ('output2', DummyType)])
            )

        self.assertEqual(signature.name, 'Bogus human-readable name')
        self.assertEqual(signature.input_artifacts,
                         {'input1': DummyType, 'input2': DummyType})
        self.assertEqual(signature.input_parameters,
                         {'param1': Int, 'param2': Int})
        self.assertEqual(signature.output_artifacts,
                         collections.OrderedDict([('output1', DummyType),
                                                  ('output2', DummyType)]))

    def test_call(self):
        signature = Signature(
            'Bogus human-readable name',
            inputs={'input1': DummyType,
                    'input2': DummyType,
                    'param1': Int,
                    'param2': Int},
            outputs=collections.OrderedDict([('output1', DummyType),
                                             ('output2', DummyType)])
            )
        self.assertEqual(signature({}, {}),
                         collections.OrderedDict([('output1', DummyType),
                                                  ('output2', DummyType)]))


if __name__ == "__main__":
    unittest.main()
