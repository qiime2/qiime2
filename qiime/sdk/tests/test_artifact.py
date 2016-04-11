# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tarfile
import tempfile
import unittest
import uuid

from qiime.plugin import Type
from qiime.sdk import Artifact, Provenance


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


class TestArtifact(unittest.TestCase):
    def setUp(self):
        # TODO standardize temporary directories created by QIIME
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-temp-')

        self.dummy_provenance = Provenance(
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

    def tearDown(self):
        self.test_dir.cleanup()

    def test_save(self):
        fp = os.path.join(self.test_dir.name, 'artifact.qtf')
        Artifact.save([-1, 42, 0, 43, 43], DummyType, self.dummy_provenance,
                      fp)

        with tarfile.open(fp, mode='r') as tar:
            fps = set(tar.getnames())
            expected = {'artifact', 'artifact/metadata.yaml',
                        'artifact/README.md', 'artifact/data',
                        'artifact/data/data.txt'}
            self.assertEqual(fps, expected)

    def test_constructor(self):
        fp = os.path.join(self.test_dir.name, 'artifact.qtf')
        Artifact.save([-1, 42, 0, 43, 43], DummyType, self.dummy_provenance,
                      fp)

        artifact = Artifact(fp)

        self.assertEqual(artifact.data, [-1, 42, 0, 43, 43])
        self.assertEqual(artifact.type, DummyType)
        self.assertEqual(artifact.provenance, self.dummy_provenance)
        self.assertIsInstance(artifact.uuid, uuid.UUID)
        self.assertEqual(artifact.uuid.version, 4)

    def test_constructor_no_provenance(self):
        fp = os.path.join(self.test_dir.name, 'artifact.qtf')
        Artifact.save([-1, 42, 0, 43, 43], DummyType, None, fp)

        artifact = Artifact(fp)

        self.assertEqual(artifact.data, [-1, 42, 0, 43, 43])
        self.assertEqual(artifact.type, DummyType)
        self.assertEqual(artifact.provenance, None)
        self.assertIsInstance(artifact.uuid, uuid.UUID)
        self.assertEqual(artifact.uuid.version, 4)


if __name__ == '__main__':
    unittest.main()
