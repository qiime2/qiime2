# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest
import uuid
import warnings

import qiime.core.type
from qiime.sdk import Result, Artifact, Visualization, Provenance
from qiime.sdk.result import ResultMetadata

from qiime.core.testing.type import FourInts
from qiime.core.testing.util import get_dummy_plugin
from qiime.core.testing.visualizer import mapping_viz


class TestResult(unittest.TestCase):
    def setUp(self):
        # Ignore the returned dummy plugin object, just run this to verify the
        # plugin exists as the tests rely on it being loaded.
        get_dummy_plugin()

        # TODO standardize temporary directories created by QIIME
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')

        self.data_dir = os.path.join(self.test_dir.name, 'viz-output')
        os.mkdir(self.data_dir)
        mapping_viz(self.data_dir,
                    {'abc': 'foo', 'def': 'bar'},
                    {'ghi': 'baz', 'jkl': 'bazz'},
                    key_label='Key', value_label='Value')

        self.provenance = Provenance(
            execution_uuid=uuid.UUID('7e909a23-21e2-44c2-be17-0723fae91dc8'),
            executor_reference=(
                'dummy_action_id. Details on plugin, version, website, etc. '
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

    def tearDown(self):
        self.test_dir.cleanup()

    def test_private_constructor(self):
        with self.assertRaisesRegex(
                NotImplementedError,
                'Result constructor.*private.*Result.load'):
            Result()

    def test_load_visualization_as_artifact(self):
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)
        visualization.save(fp)

        with self.assertRaisesRegex(
                TypeError, 'Visualization.*Artifact.load.*Visualization.load'):
            Artifact.load(fp)

    def test_load_artifact_as_visualization(self):
        artifact = Artifact._from_view(FourInts, [0, 0, 42, 1000],
                                       list, self.provenance)
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        artifact.save(fp)

        with self.assertRaisesRegex(
                TypeError, 'Artifact.*Visualization.load.*Artifact.load'):
            Visualization.load(fp)

    def test_extract_visualization_as_artifact(self):
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)
        visualization.save(fp)

        with self.assertRaisesRegex(
                TypeError, 'Artifact does not support.*Visualization'):
            Artifact.extract(fp, self.test_dir)

    def test_extract_artifact_as_visualization(self):
        artifact = Artifact._from_view(FourInts, [0, 0, 42, 1000],
                                       list, self.provenance)
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        artifact.save(fp)

        with self.assertRaisesRegex(
                TypeError, 'Visualization does not support.*FourInts'):
            Visualization.extract(fp, self.test_dir)

    def test_peek_visualization_as_artifact(self):
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)
        visualization.save(fp)

        with self.assertRaisesRegex(
                TypeError, 'Artifact does not support.*Visualization'):
            Artifact.peek(fp)

    def test_peek_artifact_as_visualization(self):
        artifact = Artifact._from_view(FourInts, [0, 0, 42, 1000],
                                       list, self.provenance)
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        artifact.save(fp)

        with self.assertRaisesRegex(
                TypeError, 'Visualization does not support.*FourInts'):
            Visualization.peek(fp)

    def test_load_artifact(self):
        saved_artifact = Artifact._from_view(FourInts, [-1, 42, 0, 43],
                                             list, self.provenance)
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        saved_artifact.save(fp)

        artifact = Result.load(fp)

        self.assertIsInstance(artifact, Artifact)
        self.assertEqual(artifact.type, FourInts)
        self.assertEqual(artifact.provenance, self.provenance)
        self.assertEqual(artifact.uuid, saved_artifact.uuid)
        self.assertEqual(artifact.view(list), [-1, 42, 0, 43])

    def test_load_visualization(self):
        saved_visualization = Visualization._from_data_dir(self.data_dir,
                                                           self.provenance)
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        saved_visualization.save(fp)

        visualization = Result.load(fp)

        self.assertIsInstance(visualization, Visualization)
        self.assertEqual(visualization.type, qiime.core.type.Visualization)
        self.assertEqual(visualization.provenance, self.provenance)
        self.assertEqual(visualization.uuid, saved_visualization.uuid)

    def test_extract_artifact(self):
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        artifact = Artifact._from_view(FourInts, [-1, 42, 0, 43],
                                       list, self.provenance)
        artifact.save(fp)

        output_dir = os.path.join(self.test_dir.name, 'artifact-extract-test')
        result_dir = Result.extract(fp, output_dir=output_dir)
        self.assertEqual(result_dir, output_dir)

        contents = [
            'artifact/VERSION',
            'artifact/metadata.yaml',
            'artifact/README.md',
            'artifact/data/file1.txt',
            'artifact/data/file2.txt',
            'artifact/data/nested/file3.txt',
            'artifact/data/nested/file4.txt']
        for fp in contents:
            expected_fp = os.path.join(output_dir, fp)
            self.assertTrue(os.path.exists(expected_fp),
                            'File %s was not extracted.' % fp)

    def test_extract_visualization(self):
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)
        visualization.save(fp)

        output_dir = os.path.join(self.test_dir.name, 'viz-extract-test')
        result_dir = Result.extract(fp, output_dir=output_dir)
        self.assertEqual(result_dir, output_dir)

        contents = [
            'visualization/VERSION',
            'visualization/metadata.yaml',
            'visualization/README.md',
            'visualization/data/index.html',
            'visualization/data/css/style.css']
        for fp in contents:
            expected_fp = os.path.join(output_dir, fp)
            self.assertTrue(os.path.exists(expected_fp),
                            'File %s was not extracted.' % fp)

    def test_peek_artifact(self):
        artifact = Artifact._from_view(FourInts, [0, 0, 42, 1000],
                                       list, self.provenance)
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        artifact.save(fp)

        metadata = Result.peek(fp)

        self.assertIsInstance(metadata, ResultMetadata)
        self.assertEqual(metadata.type, FourInts)
        self.assertEqual(metadata.provenance, self.provenance)
        self.assertEqual(metadata.uuid, artifact.uuid)

    def test_peek_visualization(self):
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization.save(fp)

        metadata = Result.peek(fp)

        self.assertIsInstance(metadata, ResultMetadata)
        self.assertEqual(metadata.type, qiime.core.type.Visualization)
        self.assertEqual(metadata.provenance, self.provenance)
        self.assertEqual(metadata.uuid, visualization.uuid)

    def test_save_artifact_warning(self):
        artifact = Artifact._from_view(FourInts, [0, 0, 42, 1000],
                                       list, self.provenance)
        fp = os.path.join(self.test_dir.name, 'artifact.zip')

        with warnings.catch_warnings(record=True) as w:
            artifact.save(fp)

            self.assertEqual(len(w), 1)
            self.assertIsInstance(w[0].message, UserWarning)
            self.assertIn(Artifact.extension, str(w[0].message))

    def test_save_visualization_warning(self):
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)
        fp = os.path.join(self.test_dir.name, 'visualization.zip')

        with warnings.catch_warnings(record=True) as w:
            visualization.save(fp)

            self.assertEqual(len(w), 1)
            self.assertIsInstance(w[0].message, UserWarning)
            self.assertIn(Visualization.extension, str(w[0].message))


if __name__ == '__main__':
    unittest.main()
