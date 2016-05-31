# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import zipfile
import tempfile
import unittest

import qiime.core.testing
from qiime.sdk import Artifact, Provenance


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

        self.artifact_with_provenance = Artifact._from_view(
            [-1, 42, 0, 43, 43], qiime.core.testing.TestType,
            self.dummy_provenance)

        self.artifact_without_provenance = Artifact._from_view(
            [-1, 42, 0, 43, 43], qiime.core.testing.TestType, None)

    def tearDown(self):
        self.test_dir.cleanup()

    def test_save(self):
        fp = os.path.join(self.test_dir.name, 'artifact.qzf')

        self.artifact_with_provenance.save(fp)

        with zipfile.ZipFile(fp, mode='r') as zf:
            fps = set(zf.namelist())
            expected = {'artifact/VERSION', 'artifact/metadata.yaml',
                        'artifact/README.md', 'artifact/data/data.txt'}
            self.assertEqual(fps, expected)

    def test_load_with_provenance(self):
        fp = os.path.join(self.test_dir.name, 'artifact.qzf')
        self.artifact_with_provenance.save(fp)

        artifact = Artifact.load(fp)

        self.assertEqual(artifact.type, qiime.core.testing.TestType)
        self.assertEqual(artifact.provenance, self.dummy_provenance)
        self.assertEqual(artifact.uuid, self.artifact_with_provenance.uuid)
        self.assertEqual(artifact.view(list), [-1, 42, 0, 43, 43])

    def test_load_without_provenance(self):
        fp = os.path.join(self.test_dir.name, 'artifact.qzf')
        self.artifact_without_provenance.save(fp)

        artifact = Artifact.load(fp)

        self.assertEqual(artifact.type, qiime.core.testing.TestType)
        self.assertEqual(artifact.provenance, None)
        self.assertEqual(artifact.uuid, self.artifact_without_provenance.uuid)
        self.assertEqual(artifact.view(list), [-1, 42, 0, 43, 43])


# TODO executing these tests via `python qiime/sdk/tests/test_artifact.py`
# for me (Jai) raises a ton of DeprecationWarnings from third-party packages I
# have installed in my conda environment (IPython, pandas, matplotlib). I
# haven't figured out if the code is doing something wrong or if I broke my
# environment somehow during testing. For now, running as
# `nosetests qiime/sdk/tests/test_artifact.py` doesn't have this issue.
if __name__ == '__main__':
    unittest.main()
