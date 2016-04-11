# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import tempfile
import collections
import concurrent
import os
import subprocess

from qiime.sdk import Workflow, Artifact, SubprocessExecutor
from qiime.plugin import Type, Int


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


def dummy_function(input1, input2, param1, param2):
    return input1 + input2 + [param1] + [param2]


class TestSubprocessExecutor(unittest.TestCase):

    def setUp(self):
        # TODO standardize temporary directories created by QIIME
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-temp-')

        self.workflow = Workflow.from_function(
            dummy_function,
            inputs={
                'input1': DummyType,
                'input2': DummyType,
                'param1': Int,
                'param2': Int
            },
            outputs=collections.OrderedDict([
                ('concatenated_inputs', DummyType)
            ]),
            name='Concatenate things',
            doc="Let's concatenate some things!"
        )

        self.artifact_fp1 = os.path.join(self.test_dir.name, 'artifact1.qtf')
        self.artifact_fp2 = os.path.join(self.test_dir.name, 'artifact2.qtf')
        Artifact.save([-1, 42, 0, 43, 43], DummyType, None, self.artifact_fp1)
        Artifact.save([1, 2, 100], DummyType, None, self.artifact_fp2)
        self.artifact_fp3 = os.path.join(self.test_dir.name, 'artifact3.qtf')

    def tearDown(self):
        self.test_dir.cleanup()

    def test_call(self):
        executor = SubprocessExecutor()
        future = executor(self.workflow,
                          input_artifact_filepaths={
                            'input1': self.artifact_fp1,
                            'input2': self.artifact_fp2
                          },
                          parameter_references={
                            'param1': 99,
                            'param2': -999,
                          },
                          output_artifact_filepaths={
                            'concatenated_inputs': self.artifact_fp3
                          })
        self.assertIsInstance(future, concurrent.futures.Future)

        completed_process = future.result()
        self.assertIsInstance(completed_process, subprocess.CompletedProcess)
        self.assertEqual(completed_process.returncode, 0)
        self.assertEqual(completed_process.stdout, b'')
        self.assertEqual(completed_process.stderr, b'')

        result = Artifact(self.artifact_fp3)
        self.assertEqual(result.data,
                         [-1, 42, 0, 43, 43, 1, 2, 100, 99, -999])
        self.assertIsNotNone(result.provenance)
        self.assertEqual(result.type, DummyType)


if __name__ == "__main__":
    unittest.main()
