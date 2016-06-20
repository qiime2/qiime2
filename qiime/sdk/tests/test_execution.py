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

import qiime.core.testing
import qiime.plugin
from qiime.sdk import Workflow, Artifact, SubprocessExecutor


def dummy_function(input1: list, input2: list,
                   param1: int, param2: int) -> list:
    return input1 + input2 + [param1] + [param2]


class TestSubprocessExecutor(unittest.TestCase):
    def setUp(self):
        # TODO standardize temporary directories created by QIIME
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-temp-')

        self.workflow = Workflow.from_function(
            dummy_function,
            inputs={
                'input1': qiime.core.testing.TestType,
                'input2': qiime.core.testing.TestType
            },
            parameters={
                'param1': qiime.plugin.Int,
                'param2': qiime.plugin.Int
            },
            outputs=collections.OrderedDict([
                ('concatenated_inputs', qiime.core.testing.TestType)
            ]),
            name='Concatenate things',
            doc="Let's concatenate some things!"
        )

        self.artifact_fp1 = os.path.join(self.test_dir.name, 'artifact1.qza')
        self.artifact_fp2 = os.path.join(self.test_dir.name, 'artifact2.qza')
        self.artifact_fp3 = os.path.join(self.test_dir.name, 'artifact3.qza')

        artifact = Artifact._from_view(
            [-1, 42, 0, 43, 43], qiime.core.testing.TestType, None)
        artifact.save(self.artifact_fp1)

        artifact = Artifact._from_view(
            [1, 2, 100], qiime.core.testing.TestType, None)
        artifact.save(self.artifact_fp2)

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

        artifact = Artifact.load(self.artifact_fp3)

        self.assertIsNotNone(artifact.provenance)
        self.assertEqual(artifact.type, qiime.core.testing.TestType)
        self.assertEqual(artifact.view(list),
                         [-1, 42, 0, 43, 43, 1, 2, 100, 99, -999])


# TODO executing these tests via `python qiime/sdk/tests/test_execution.py`
# doesn't work, currently must be run as
# `nosetests qiime/sdk/tests/test_execution.py`
if __name__ == "__main__":
    unittest.main()
