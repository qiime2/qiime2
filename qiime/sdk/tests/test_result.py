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

import qiime.core.type
from qiime.sdk import Result, Artifact, Visualization
from qiime.sdk.result import ResultMetadata
import qiime.core.archive as archive

from qiime.core.testing.type import FourInts
from qiime.core.testing.util import get_dummy_plugin, ArchiveTestingMixin
from qiime.core.testing.visualizer import mapping_viz


class TestResult(unittest.TestCase, ArchiveTestingMixin):
    def make_provenance_capture(self):
        # You can't actually import a visualization, but I won't tell
        # visualization if you don't...
        return archive.ImportProvenanceCapture()

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

    def tearDown(self):
        self.test_dir.cleanup()

    def test_private_constructor(self):
        with self.assertRaisesRegex(
                NotImplementedError,
                'Result constructor.*private.*Result.load'):
            Result()

    def test_load_artifact(self):
        saved_artifact = Artifact.import_data(FourInts, [-1, 42, 0, 43])
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        saved_artifact.save(fp)

        artifact = Result.load(fp)

        self.assertIsInstance(artifact, Artifact)
        self.assertEqual(artifact.type, FourInts)
        self.assertEqual(artifact.uuid, saved_artifact.uuid)
        self.assertEqual(artifact.view(list), [-1, 42, 0, 43])

    def test_load_visualization(self):
        saved_visualization = Visualization._from_data_dir(
             self.data_dir, self.make_provenance_capture())
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        saved_visualization.save(fp)

        visualization = Result.load(fp)

        self.assertIsInstance(visualization, Visualization)
        self.assertEqual(visualization.type, qiime.core.type.Visualization)
        self.assertEqual(visualization.uuid, saved_visualization.uuid)

    def test_extract_artifact(self):
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        artifact = Artifact.import_data(FourInts, [-1, 42, 0, 43])
        artifact.save(fp)

        root_dir = str(artifact.uuid)
        output_dir = os.path.join(self.test_dir.name, 'artifact-extract-test')
        result_dir = Result.extract(fp, output_dir=output_dir)
        self.assertEqual(result_dir, os.path.join(output_dir, root_dir))

        expected = {
            'VERSION',
            'metadata.yaml',
            'data/file1.txt',
            'data/file2.txt',
            'data/nested/file3.txt',
            'data/nested/file4.txt',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/action/action.yaml'
        }

        self.assertExtractedArchiveMembers(output_dir, root_dir, expected)

    def test_extract_visualization(self):
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(
             self.data_dir, self.make_provenance_capture())
        visualization.save(fp)

        root_dir = str(visualization.uuid)
        output_dir = os.path.join(self.test_dir.name, 'viz-extract-test')
        result_dir = Result.extract(fp, output_dir=output_dir)
        self.assertEqual(result_dir, os.path.join(output_dir, root_dir))

        expected = {
            'VERSION',
            'metadata.yaml',
            'data/index.html',
            'data/css/style.css',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/action/action.yaml'
        }

        self.assertExtractedArchiveMembers(output_dir, root_dir, expected)

    def test_peek_artifact(self):
        artifact = Artifact.import_data(FourInts, [0, 0, 42, 1000])
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        artifact.save(fp)

        metadata = Result.peek(fp)

        self.assertIsInstance(metadata, ResultMetadata)
        self.assertEqual(metadata.type, 'FourInts')
        self.assertEqual(metadata.uuid, str(artifact.uuid))
        self.assertEqual(metadata.format, 'FourIntsDirectoryFormat')

    def test_peek_visualization(self):
        visualization = Visualization._from_data_dir(
             self.data_dir, self.make_provenance_capture())
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization.save(fp)

        metadata = Result.peek(fp)

        self.assertIsInstance(metadata, ResultMetadata)
        self.assertEqual(metadata.type, 'Visualization')
        self.assertEqual(metadata.uuid, str(visualization.uuid))
        self.assertIsNone(metadata.format)

    def test_save_artifact_auto_extension(self):
        artifact = Artifact.import_data(FourInts, [0, 0, 42, 1000])

        # No extension.
        fp = os.path.join(self.test_dir.name, 'artifact')
        obs_fp = artifact.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.qza')

        # Wrong extension.
        fp = os.path.join(self.test_dir.name, 'artifact.zip')
        obs_fp = artifact.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.zip.qza')

        # Correct extension.
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        obs_fp = artifact.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.qza')

    def test_save_visualization_auto_extension(self):
        visualization = Visualization._from_data_dir(
             self.data_dir, self.make_provenance_capture())

        # No extension.
        fp = os.path.join(self.test_dir.name, 'visualization')
        obs_fp = visualization.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.qzv')

        # Wrong extension.
        fp = os.path.join(self.test_dir.name, 'visualization.zip')
        obs_fp = visualization.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.zip.qzv')

        # Correct extension.
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        obs_fp = visualization.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.qzv')


if __name__ == '__main__':
    unittest.main()
