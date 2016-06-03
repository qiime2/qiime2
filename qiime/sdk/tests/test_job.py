# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import uuid

from qiime.sdk.job import Job


class TestJob(unittest.TestCase):

    def test_constructor(self):
        job = Job(code="import this",
                  uuid=uuid.UUID('7e909a23-21e2-44c2-be17-0723fae91dc8'),
                  input_artifact_filepaths={'input1': 'input1.qzf'},
                  parameter_references={'param1': 42},
                  output_artifact_filepaths={'output1': 'output1.qzf'})

        self.assertEqual(job.code, "import this")
        self.assertEqual(job.uuid,
                         uuid.UUID('7e909a23-21e2-44c2-be17-0723fae91dc8'))
        self.assertEqual(job.input_artifact_filepaths,
                         {'input1': 'input1.qzf'})
        self.assertEqual(job.parameter_references, {'param1': 42})
        self.assertEqual(job.output_artifact_filepaths,
                         {'output1': 'output1.qzf'})


if __name__ == "__main__":
    unittest.main()
