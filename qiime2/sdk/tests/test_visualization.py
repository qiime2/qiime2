# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest
import uuid
import collections
import pathlib

import qiime2.core.type
from qiime2.sdk import Visualization
from qiime2.sdk.result import ResultMetadata
import qiime2.core.archive as archive

from qiime2.core.testing.visualizer import (
    mapping_viz, most_common_viz, multi_html_viz)
from qiime2.core.testing.util import ArchiveTestingMixin


class TestVisualization(unittest.TestCase, ArchiveTestingMixin):
    def make_provenance_capture(self):
        # You can't actually import a visualization, but I won't tell
        # visualization if you don't...
        return archive.ImportProvenanceCapture()

    def setUp(self):
        # TODO standardize temporary directories created by QIIME 2
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')

        # Using `mapping_viz` because it produces multiple files, including a
        # nested directory.
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
                'Visualization constructor.*private.*Visualization.load'):
            Visualization()

    # Note on testing strategy below: many of the tests for `_from_data_dir`
    # and `load` are similar, with the exception that when `load`ing, the
    # visualization's UUID is known so more specific assertions can be
    # performed. While these tests appear somewhat redundant, they are
    # important because they exercise the same operations on Visualization
    # objects constructed from different sources, whose codepaths have very
    # different internal behavior. This internal behavior could be tested
    # explicitly but it is safer to test the public API behavior (e.g. as a
    # user would interact with the object) in case the internals change.

    def test_from_data_dir(self):
        visualization = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())

        self.assertEqual(visualization.type, qiime2.core.type.Visualization)
        self.assertIsInstance(visualization.uuid, uuid.UUID)

    def test_from_data_dir_and_save(self):
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())

        visualization.save(fp)

        root_dir = str(visualization.uuid)
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

        self.assertArchiveMembers(fp, root_dir, expected)

    def test_load(self):
        saved_visualization = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        saved_visualization.save(fp)

        visualization = Visualization.load(fp)

        self.assertEqual(visualization.type, qiime2.core.type.Visualization)
        self.assertEqual(visualization.uuid, saved_visualization.uuid)

    def test_load_and_save(self):
        fp1 = os.path.join(self.test_dir.name, 'visualization1.qzv')
        fp2 = os.path.join(self.test_dir.name, 'visualization2.qzv')
        visualization = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())
        visualization.save(fp1)

        visualization = Visualization.load(fp1)
        # Overwriting its source file works.
        visualization.save(fp1)
        # Saving to a new file works.
        visualization.save(fp2)

        root_dir = str(visualization.uuid)
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

        self.assertArchiveMembers(fp1, root_dir, expected)

        root_dir = str(visualization.uuid)
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

        self.assertArchiveMembers(fp2, root_dir, expected)

    def test_roundtrip(self):
        fp1 = os.path.join(self.test_dir.name, 'visualization1.qzv')
        fp2 = os.path.join(self.test_dir.name, 'visualization2.qzv')
        visualization = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())
        visualization.save(fp1)

        visualization1 = Visualization.load(fp1)
        visualization1.save(fp2)
        visualization2 = Visualization.load(fp2)

        self.assertEqual(visualization1.type, visualization2.type)
        self.assertEqual(visualization1.uuid, visualization2.uuid)

    def test_load_with_archive_filepath_modified(self):
        # Save a visualization for use in the following test case.
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        Visualization._from_data_dir(self.data_dir,
                                     self.make_provenance_capture()).save(fp)

        # Load the visualization from a filepath then save a different
        # visualization to the same filepath. Assert that both visualizations
        # access the correct data.
        #
        # `load` used to be lazy, only extracting data when it needed to (e.g.
        # when `save` or `get_index_paths` was called). This was buggy as the
        # filepath could have been deleted, or worse, modified to contain a
        # different .qzv file. Thus, the wrong archive could be extracted on
        # demand, or the archive could be missing altogether. There isn't an
        # easy cross-platform compatible way to solve this problem, so
        # Visualization.load is no longer lazy and always extracts its data
        # immediately. The real motivation for lazy loading was for quick
        # inspection of archives without extracting/copying data, so that API
        # is now provided through Visualization.peek.
        visualization1 = Visualization.load(fp)

        new_data_dir = os.path.join(self.test_dir.name, 'viz-output2')
        os.mkdir(new_data_dir)
        most_common_viz(new_data_dir, collections.Counter(range(42)))

        Visualization._from_data_dir(new_data_dir,
                                     self.make_provenance_capture()).save(fp)
        visualization2 = Visualization.load(fp)

        self.assertEqual(visualization1.get_index_paths(),
                         {'html': 'data/index.html'})
        self.assertEqual(visualization2.get_index_paths(),
                         {'html': 'data/index.html', 'tsv': 'data/index.tsv'})

    def test_extract(self):
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())
        visualization.save(fp)

        root_dir = str(visualization.uuid)
        # pathlib normalizes away the `.`, it doesn't matter, but this is the
        # implementation we're using, so let's test against that assumption.
        output_dir = pathlib.Path(self.test_dir.name) / 'viz-extract-test'
        result_dir = Visualization.extract(fp, output_dir=output_dir)
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

    def test_get_index_paths_single_load(self):
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())
        visualization.save(fp)
        visualization = Visualization.load(fp)

        actual = visualization.get_index_paths()
        expected = {'html': 'data/index.html'}
        self.assertEqual(actual, expected)

    def test_get_index_paths_single_from_data_dir(self):
        visualization = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())

        actual = visualization.get_index_paths()
        expected = {'html': 'data/index.html'}
        self.assertEqual(actual, expected)

    def test_get_index_paths_multiple_load(self):
        data_dir = os.path.join(self.test_dir.name, 'mc-viz-output1')
        os.mkdir(data_dir)
        most_common_viz(data_dir,
                        collections.Counter(range(42)))
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(
            data_dir, self.make_provenance_capture())
        visualization.save(fp)
        visualization = Visualization.load(fp)

        actual = visualization.get_index_paths()
        expected = {'html': 'data/index.html',
                    'tsv': 'data/index.tsv'}
        self.assertEqual(actual, expected)

    def test_get_index_paths_multiple_from_data_dir(self):
        data_dir = os.path.join(self.test_dir.name, 'mc-viz-output2')
        os.mkdir(data_dir)
        most_common_viz(data_dir, collections.Counter(range(42)))
        visualization = Visualization._from_data_dir(
            data_dir, self.make_provenance_capture())

        actual = visualization.get_index_paths()
        expected = {'html': 'data/index.html',
                    'tsv': 'data/index.tsv'}
        self.assertEqual(actual, expected)

    def test_get_index_paths_multiple_html_load(self):
        data_dir = os.path.join(self.test_dir.name, 'multi-html-viz1')
        os.mkdir(data_dir)
        multi_html_viz(data_dir, [1, 42])

        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(
            data_dir, self.make_provenance_capture())
        visualization.save(fp)
        visualization = Visualization.load(fp)

        with self.assertRaises(ValueError):
            visualization.get_index_paths()

    def test_get_index_paths_multiple_html_from_data_dir(self):
        data_dir = os.path.join(self.test_dir.name, 'multi-html-viz2')
        os.mkdir(data_dir)
        multi_html_viz(data_dir, [1, 42])

        visualization = Visualization._from_data_dir(
            data_dir, self.make_provenance_capture())

        with self.assertRaises(ValueError):
            visualization.get_index_paths()

    def test_get_index_paths_relative_false(self):
        data_dir = os.path.join(self.test_dir.name, 'mc-viz-output2')
        os.mkdir(data_dir)
        most_common_viz(data_dir, collections.Counter(range(42)))
        visualization = Visualization._from_data_dir(
            data_dir, self.make_provenance_capture())

        def get_abs_path(rel):
            return str(visualization._archiver.root_dir / rel)
        actual = visualization.get_index_paths(relative=False)
        expected = {'html': get_abs_path('data/index.html'),
                    'tsv': get_abs_path('data/index.tsv')}
        self.assertEqual(actual, expected)

    def test_peek(self):
        visualization = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization.save(fp)

        metadata = Visualization.peek(fp)

        self.assertIsInstance(metadata, ResultMetadata)
        self.assertEqual(metadata.type, 'Visualization')
        self.assertEqual(metadata.uuid, str(visualization.uuid))
        self.assertIsNone(metadata.format)

    def test_eq_identity(self):
        visualization = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())

        self.assertEqual(visualization, visualization)

    def test_eq_same_uuid(self):
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization1 = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())
        visualization1.save(fp)

        visualization2 = Visualization.load(fp)

        self.assertEqual(visualization1, visualization2)

    def test_ne_same_data_different_uuid(self):
        visualization1 = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())
        visualization2 = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())

        self.assertNotEqual(visualization1, visualization2)

    def test_ne_different_data_different_uuid(self):
        visualization1 = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())

        data_dir = os.path.join(self.test_dir.name, 'mc-viz-output1')
        os.mkdir(data_dir)
        most_common_viz(data_dir,
                        collections.Counter(range(42)))
        visualization2 = Visualization._from_data_dir(
            data_dir, self.make_provenance_capture())

        self.assertNotEqual(visualization1, visualization2)

    def test_ne_subclass_same_uuid(self):
        class VisualizationSubclass(Visualization):
            pass

        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization1 = VisualizationSubclass._from_data_dir(
            self.data_dir, self.make_provenance_capture())
        visualization1.save(fp)

        visualization2 = Visualization.load(fp)

        self.assertNotEqual(visualization1, visualization2)
        self.assertNotEqual(visualization2, visualization1)

    def test_ne_different_type_same_uuid(self):
        visualization = Visualization._from_data_dir(
            self.data_dir, self.make_provenance_capture())

        class Faker:
            @property
            def uuid(self):
                return visualization.uuid

        faker = Faker()

        self.assertNotEqual(visualization, faker)


if __name__ == '__main__':
    unittest.main()
