# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest
import pathlib

import qiime2.core.type
from qiime2.sdk import Result, Artifact, Visualization, ResultCollection
from qiime2.sdk.result import ResultMetadata
import qiime2.core.archive as archive
import qiime2.core.exceptions as exceptions

from qiime2.core.testing.type import (FourInts, SingleInt, IntSequence1,
                                      IntSequence2)
from qiime2.core.testing.util import get_dummy_plugin, ArchiveTestingMixin
from qiime2.core.testing.visualizer import mapping_viz
from qiime2.core.util import set_permissions, OTHER_NO_WRITE


class TestResult(unittest.TestCase, ArchiveTestingMixin):
    def make_provenance_capture(self):
        # You can't actually import a visualization, but I won't tell
        # visualization if you don't...
        return archive.ImportProvenanceCapture()

    def setUp(self):
        # Ignore the returned dummy plugin object, just run this to verify the
        # plugin exists as the tests rely on it being loaded.
        get_dummy_plugin()

        # TODO standardize temporary directories created by QIIME 2
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
        self.assertEqual(visualization.type, qiime2.core.type.Visualization)
        self.assertEqual(visualization.uuid, saved_visualization.uuid)

    def test_extract_artifact(self):
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        artifact = Artifact.import_data(FourInts, [-1, 42, 0, 43])
        artifact.save(fp)

        root_dir = str(artifact.uuid)
        # pathlib normalizes away the `.`, it doesn't matter, but this is the
        # implementation we're using, so let's test against that assumption.
        output_dir = pathlib.Path(self.test_dir.name) / 'artifact-extract-test'
        result_dir = Result.extract(fp, output_dir=output_dir)
        self.assertEqual(result_dir, str(output_dir / root_dir))

        expected = {
            'VERSION',
            'checksums.md5',
            'metadata.yaml',
            'data/file1.txt',
            'data/file2.txt',
            'data/nested/file3.txt',
            'data/nested/file4.txt',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/citations.bib',
            'provenance/action/action.yaml'
        }

        self.assertExtractedArchiveMembers(output_dir, root_dir, expected)

    def test_extract_visualization(self):
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(
             self.data_dir, self.make_provenance_capture())
        visualization.save(fp)

        root_dir = str(visualization.uuid)
        output_dir = pathlib.Path(self.test_dir.name) / 'viz-extract-test'
        result_dir = Result.extract(fp, output_dir=output_dir)
        self.assertEqual(result_dir, str(output_dir / root_dir))

        expected = {
            'VERSION',
            'checksums.md5',
            'metadata.yaml',
            'data/index.html',
            'data/css/style.css',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/citations.bib',
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

        # Filename & extension endswith is matching (default).
        fp = os.path.join(self.test_dir.name, 'artifactqza')
        obs_fp = artifact.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifactqza.qza')

        # Filename & extension endswith is matching (non-default).
        fp = os.path.join(self.test_dir.name, 'artifacttxt')
        obs_fp = artifact.save(fp, 'txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifacttxt.txt')

        # No period in filename; no period in extension.
        fp = os.path.join(self.test_dir.name, 'artifact')
        obs_fp = artifact.save(fp, 'txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.txt')

        # No period in filename; multiple periods in extension.
        fp = os.path.join(self.test_dir.name, 'artifact')
        obs_fp = artifact.save(fp, '..txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.txt')

        # Single period in filename; no period in extension.
        fp = os.path.join(self.test_dir.name, 'artifact.')
        obs_fp = artifact.save(fp, 'txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.txt')

        # Single period in filename; single period in extension.
        fp = os.path.join(self.test_dir.name, 'artifact.')
        obs_fp = artifact.save(fp, '.txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.txt')

        # Single period in filename; multiple periods in extension.
        fp = os.path.join(self.test_dir.name, 'artifact.')
        obs_fp = artifact.save(fp, '..txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.txt')

        # Multiple periods in filename; single period in extension.
        fp = os.path.join(self.test_dir.name, 'artifact..')
        obs_fp = artifact.save(fp, '.txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.txt')

        # Multiple periods in filename; multiple periods in extension.
        fp = os.path.join(self.test_dir.name, 'artifact..')
        obs_fp = artifact.save(fp, '..txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.txt')

        # No extension in filename; no extension input.
        fp = os.path.join(self.test_dir.name, 'artifact')
        obs_fp = artifact.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.qza')

        # No extension in filename; different extension input.
        fp = os.path.join(self.test_dir.name, 'artifact')
        obs_fp = artifact.save(fp, '.txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.txt')

        # No extension in filename; default extension input.
        fp = os.path.join(self.test_dir.name, 'artifact')
        obs_fp = artifact.save(fp, '.qza')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.qza')

        # Different extension in filename; no extension input.
        fp = os.path.join(self.test_dir.name, 'artifact.zip')
        obs_fp = artifact.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.zip.qza')

        # Different extension in filename;
        # Different extension input (non-matching).
        fp = os.path.join(self.test_dir.name, 'artifact.zip')
        obs_fp = artifact.save(fp, '.txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.zip.txt')

        # Different extension in filename;
        # Different extension input (matching).
        fp = os.path.join(self.test_dir.name, 'artifact.zip')
        obs_fp = artifact.save(fp, '.zip')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.zip')

        # Different extension in filename; default extension input.
        fp = os.path.join(self.test_dir.name, 'artifact.zip')
        obs_fp = artifact.save(fp, '.qza')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.zip.qza')

        # Default extension in filename; no extension input.
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        obs_fp = artifact.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.qza')

        # Default extension in filename; different extension input.
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        obs_fp = artifact.save(fp, '.txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.qza.txt')

        # Default extension in filename; default extension input.
        fp = os.path.join(self.test_dir.name, 'artifact.qza')
        obs_fp = artifact.save(fp, '.qza')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'artifact.qza')

    def test_save_visualization_auto_extension(self):
        visualization = Visualization._from_data_dir(
             self.data_dir, self.make_provenance_capture())

        # Filename & extension endswith is matching (default).
        fp = os.path.join(self.test_dir.name, 'visualizationqzv')
        obs_fp = visualization.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualizationqzv.qzv')

        # Filename & extension endswith is matching (non-default).
        fp = os.path.join(self.test_dir.name, 'visualizationtxt')
        obs_fp = visualization.save(fp, 'txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualizationtxt.txt')

        # No period in filename; no period in extension.
        fp = os.path.join(self.test_dir.name, 'visualization')
        obs_fp = visualization.save(fp, 'txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.txt')

        # No period in filename; multiple periods in extension.
        fp = os.path.join(self.test_dir.name, 'visualization')
        obs_fp = visualization.save(fp, '..txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.txt')

        # Single period in filename; no period in extension.
        fp = os.path.join(self.test_dir.name, 'visualization.')
        obs_fp = visualization.save(fp, 'txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.txt')

        # Single period in filename; single period in extension.
        fp = os.path.join(self.test_dir.name, 'visualization.')
        obs_fp = visualization.save(fp, '.txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.txt')

        # Single period in filename; multiple periods in extension.
        fp = os.path.join(self.test_dir.name, 'visualization.')
        obs_fp = visualization.save(fp, '..txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.txt')

        # Multiple periods in filename; single period in extension.
        fp = os.path.join(self.test_dir.name, 'visualization..')
        obs_fp = visualization.save(fp, '.txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.txt')

        # Multiple periods in filename; multiple periods in extension.
        fp = os.path.join(self.test_dir.name, 'visualization..')
        obs_fp = visualization.save(fp, '..txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.txt')

        # No extension in filename; no extension input.
        fp = os.path.join(self.test_dir.name, 'visualization')
        obs_fp = visualization.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.qzv')

        # No extension in filename; different extension input.
        fp = os.path.join(self.test_dir.name, 'visualization')
        obs_fp = visualization.save(fp, '.txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.txt')

        # No extension in filename; default extension input.
        fp = os.path.join(self.test_dir.name, 'visualization')
        obs_fp = visualization.save(fp, '.qzv')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.qzv')

        # Different extension in filename; no extension input.
        fp = os.path.join(self.test_dir.name, 'visualization.zip')
        obs_fp = visualization.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.zip.qzv')

        # Different extension in filename;
        # Different extension input (non-matching).
        fp = os.path.join(self.test_dir.name, 'visualization.zip')
        obs_fp = visualization.save(fp, '.txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.zip.txt')

        # Different extension in filename;
        # Different extension input (matching).
        fp = os.path.join(self.test_dir.name, 'visualization.zip')
        obs_fp = visualization.save(fp, '.zip')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.zip')

        # Different extension in filename; default extension input.
        fp = os.path.join(self.test_dir.name, 'visualization.zip')
        obs_fp = visualization.save(fp, '.qzv')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.zip.qzv')

        # Default extension in filename; no extension input.
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        obs_fp = visualization.save(fp)
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.qzv')

        # Default extension in filename; different extension input.
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        obs_fp = visualization.save(fp, '.txt')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.qzv.txt')

        # Default extension in filename; default extension input.
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        obs_fp = visualization.save(fp, '.qzv')
        obs_filename = os.path.basename(obs_fp)

        self.assertEqual(obs_filename, 'visualization.qzv')

    def test_import_data_single_dirfmt_to_single_dirfmt(self):
        temp_data_dir = os.path.join(self.test_dir.name, 'import')
        os.mkdir(temp_data_dir)

        with open(os.path.join(temp_data_dir, 'ints.txt'), 'w') as fh:
            fh.write("1\n2\n3\n")

        qiime2.Artifact.import_data('IntSequence2', temp_data_dir,
                                    view_type="IntSequenceDirectoryFormat")

    def test_artifact_has_metadata_true(self):
        A = Artifact.import_data('Mapping', {'a': '1', 'b': '2'})
        self.assertTrue(A.has_metadata())

    def test_artifact_has_metadata_false(self):
        A = Artifact.import_data('IntSequence1', [1, 2, 3, 4])
        self.assertFalse(A.has_metadata())

    def test_validate_artifact_good(self):
        artifact = Artifact.import_data('IntSequence1', [1, 2, 3, 4])

        artifact.validate()
        self.assertTrue(True)  # Checkpoint

    def test_validate_artifact_bad(self):
        artifact = Artifact.import_data('IntSequence1', [1, 2, 3, 4])
        # We set everything in the artifact to be read-only. This test needs to
        # mimic if the user were to somehow write it anyway, so we set write
        # for self and group
        set_permissions(artifact._archiver.root_dir, OTHER_NO_WRITE,
                        OTHER_NO_WRITE)

        with (artifact._archiver.root_dir / 'extra.file').open('w') as fh:
            fh.write('uh oh')

        with self.assertRaisesRegex(exceptions.ValidationError,
                                    r'extra\.file'):
            artifact.validate()

    def test_validate_vizualization_good(self):
        visualization = Visualization._from_data_dir(
             self.data_dir, self.make_provenance_capture())

        visualization.validate()
        self.assertTrue(True)  # Checkpoint

    def test_validate_vizualization_bad(self):
        visualization = Visualization._from_data_dir(
             self.data_dir, self.make_provenance_capture())

        # We set everything in the artifact to be read-only. This test needs to
        # mimic if the user were to somehow write it anyway, so we set write
        # for self and group
        set_permissions(visualization._archiver.root_dir, OTHER_NO_WRITE,
                        OTHER_NO_WRITE)

        with (visualization._archiver.root_dir / 'extra.file').open('w') as fh:
            fh.write('uh oh')

        with self.assertRaisesRegex(exceptions.ValidationError,
                                    r'extra\.file'):
            visualization.validate()


class TestResultCollection(unittest.TestCase):
    def setUp(self):
        # Ignore the returned dummy plugin object, just run this to verify the
        # plugin exists as the tests rely on it being loaded.
        get_dummy_plugin()

        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        self.output_fp = os.path.join(self.test_dir.name, 'output')

        self.collection = ResultCollection(
            {'foo': Artifact.import_data(SingleInt, 0),
             'bar': Artifact.import_data(SingleInt, 1)})

    def tearDown(self):
        self.test_dir.cleanup()

    def test_roundtrip_ordered_collection(self):
        self.collection.save(self.output_fp)

        foo = Artifact.load(os.path.join(self.output_fp, 'foo.qza'))
        bar = Artifact.load(os.path.join(self.output_fp, 'bar.qza'))

        self.assertEqual(foo.view(int), 0)
        self.assertEqual(bar.view(int), 1)

        with open(os.path.join(self.output_fp, '.order')) as fh:
            self.assertEqual(fh.read(), 'foo\nbar\n')

        read_collection = ResultCollection.load(self.output_fp)
        self.assertEqual(self.collection, read_collection)

    def test_roundtrip_unordered_collection(self):
        self.collection.save(self.output_fp)
        os.remove(os.path.join(self.output_fp, '.order'))

        foo = Artifact.load(os.path.join(self.output_fp, 'foo.qza'))
        bar = Artifact.load(os.path.join(self.output_fp, 'bar.qza'))

        self.assertEqual(foo.view(int), 0)
        self.assertEqual(bar.view(int), 1)

        with self.assertWarnsRegex(
                UserWarning, f"The directory '{self.output_fp}' does not "
                "contain a .order file"):
            read_collection = ResultCollection.load(self.output_fp)

        self.assertEqual(
            set(self.collection.items()), set(read_collection.items()))

    def test_type_normal_collection(self):
        self.assertEqual(
            self.collection.type, qiime2.core.type.Collection[SingleInt])

    def test_type_weird_collection(self):
        weird_collection = ResultCollection({
            'foo': Artifact.import_data(SingleInt, 0),
            'bar': Artifact.import_data(FourInts, [1, 2, 3, 4]),
            'baz': Artifact.import_data(IntSequence1, [5, 6, 7]),
            'qux': Artifact.import_data(IntSequence2, [8, 9, 10])})

        self.assertEqual(
            weird_collection.type,
            qiime2.core.type.Collection[SingleInt | FourInts | IntSequence1 |
                                        IntSequence2])

    def test_collection_order_file_contains_nonexistent_key(self):
        BAD_KEY = 'NonexistentKey'

        self.collection.save(self.output_fp)
        order_fp = os.path.join(self.output_fp, '.order')

        with open(order_fp, 'a') as order_fh:
            order_fh.write(BAD_KEY)

        foo = Artifact.load(os.path.join(self.output_fp, 'foo.qza'))
        bar = Artifact.load(os.path.join(self.output_fp, 'bar.qza'))

        self.assertEqual(foo.view(int), 0)
        self.assertEqual(bar.view(int), 1)

        with self.assertRaisesRegex(
                ValueError, f"The Result '{BAD_KEY}' is referenced in the "
                            "order file but does not exist"):
            ResultCollection.load(self.output_fp)

    def test_collection_non_str_keys(self):
        with self.assertRaisesRegex(
                TypeError, 'All ResultCollection keys must be strings'):
            ResultCollection({1: 0})


if __name__ == '__main__':
    unittest.main()
