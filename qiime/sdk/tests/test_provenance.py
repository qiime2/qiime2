# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime.sdk import Provenance


class TestProvenance(unittest.TestCase):

    def test_constructor(self):
        provenance = Provenance(
            job_uuid='7e909a23-21e2-44c2-be17-0723fae91dc8',
            artifact_uuids={
                'input1': 'f16ca3d0-fe83-4b1e-8eea-7e35db3f6b0f',
                'input2': '908dece5-db23-4562-ad03-876bb5750145',
            },
            parameters={
                'param1': 'abc',
                'param2': 100,
            },
            workflow_reference=(
                'dummy workflow reference, see '
                'https://github.com/biocore/qiime2/issues/26'
            )
        )
        self.assertEqual(provenance.job_uuid,
                         '7e909a23-21e2-44c2-be17-0723fae91dc8')
        self.assertEqual(provenance.artifact_uuids,
                         {'input1': 'f16ca3d0-fe83-4b1e-8eea-7e35db3f6b0f',
                          'input2': '908dece5-db23-4562-ad03-876bb5750145'})
        self.assertEqual(provenance.parameters,
                         {'param1': 'abc', 'param2': 100})
        self.assertEqual(provenance.workflow_reference,
                         'dummy workflow reference, see '
                         'https://github.com/biocore/qiime2/issues/26')


if __name__ == "__main__":
    unittest.main()
