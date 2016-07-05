# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import uuid

import qiime.sdk


class TestProvenance(unittest.TestCase):
    def setUp(self):
        self.provenance = qiime.sdk.Provenance(
            execution_uuid=uuid.UUID('7e909a23-21e2-44c2-be17-0723fae91dc8'),
            executor_reference=(
                'concatenate_ints. Details on plugin, version, website, etc. '
                'will also be included, see '
                'https://github.com/biocore/qiime2/issues/26'
            ),
            artifact_uuids={
                'input1': uuid.UUID('f16ca3d0-fe83-4b1e-8eea-7e35db3f6b0f'),
                'input2': uuid.UUID('908dece5-db23-4562-ad03-876bb5750145')
            },
            parameter_references={
                'param1': 'abc',
                'param2': '100'
            }
        )

    def test_execution_uuid(self):
        self.assertEqual(self.provenance.execution_uuid,
                         uuid.UUID('7e909a23-21e2-44c2-be17-0723fae91dc8'))

    def test_executor_reference(self):
        self.assertEqual(
            self.provenance.executor_reference,
            'concatenate_ints. Details on plugin, version, website, etc. will '
            'also be included, see https://github.com/biocore/qiime2/issues/26'
        )

    def test_artifact_uuids(self):
        self.assertEqual(self.provenance.artifact_uuids, {
            'input1': uuid.UUID('f16ca3d0-fe83-4b1e-8eea-7e35db3f6b0f'),
            'input2': uuid.UUID('908dece5-db23-4562-ad03-876bb5750145')
        })

    def test_parameter_references(self):
        self.assertEqual(self.provenance.parameter_references, {
            'param1': 'abc',
            'param2': '100'
        })


if __name__ == "__main__":
    unittest.main()
