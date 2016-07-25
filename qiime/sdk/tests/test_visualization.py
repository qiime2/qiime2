# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources
import tempfile
import unittest
import uuid
import zipfile
import collections

import qiime.core.type
from qiime.sdk import Visualization, Provenance
from qiime.sdk.result import ResultMetadata

from qiime.core.testing.visualizer import (
    mapping_viz, most_common_viz, multi_html_viz)


class TestVisualization(unittest.TestCase):
    def setUp(self):
        # TODO standardize temporary directories created by QIIME
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')

        # Using `mapping_viz` because it produces multiple files, including a
        # nested directory.
        self.data_dir = os.path.join(self.test_dir.name, 'viz-output')
        os.mkdir(self.data_dir)
        mapping_viz(self.data_dir,
                    {'abc': 'foo', 'def': 'bar'},
                    {'ghi': 'baz', 'jkl': 'bazz'},
                    key_label='Key', value_label='Value')

        self.provenance = Provenance(
            execution_uuid=uuid.UUID('7e909a23-21e2-44c2-be17-0723fae91dc8'),
            executor_reference=(
                'mapping_viz. Details on plugin, version, website, etc. '
                'will also be included, see '
                'https://github.com/biocore/qiime2/issues/26'
            ),
            artifact_uuids={
                'mapping1': uuid.UUID('f16ca3d0-fe83-4b1e-8eea-7e35db3f6b0f'),
                'mapping2': uuid.UUID('908dece5-db23-4562-ad03-876bb5750145')
            },
            parameter_references={
                'key_label': 'Key',
                'value_label': 'Value'
            }
        )

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
        visualization = Visualization._from_data_dir(self.data_dir, None)

        self.assertEqual(visualization.type, qiime.core.type.Visualization)
        self.assertIsNone(visualization.provenance)
        # We don't know what the UUID is because it's generated within
        # Visualization._from_data_dir.
        self.assertIsInstance(visualization.uuid, uuid.UUID)

    def test_from_data_dir_with_provenance(self):
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)

        self.assertEqual(visualization.type, qiime.core.type.Visualization)
        self.assertEqual(visualization.provenance, self.provenance)
        self.assertIsInstance(visualization.uuid, uuid.UUID)

    def test_from_data_dir_and_save(self):
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)

        visualization.save(fp)

        with zipfile.ZipFile(fp, mode='r') as zf:
            fps = set(zf.namelist())
            expected = {
                'visualization/VERSION',
                'visualization/metadata.yaml',
                'visualization/README.md',
                'visualization/data/index.html',
                'visualization/data/css/style.css'
            }
            self.assertEqual(fps, expected)

    def test_load(self):
        saved_visualization = Visualization._from_data_dir(self.data_dir, None)
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        saved_visualization.save(fp)

        visualization = Visualization.load(fp)

        self.assertEqual(visualization.type, qiime.core.type.Visualization)
        self.assertIsNone(visualization.provenance)
        self.assertEqual(visualization.uuid, saved_visualization.uuid)

    def test_load_with_provenance(self):
        saved_visualization = Visualization._from_data_dir(self.data_dir,
                                                           self.provenance)
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        saved_visualization.save(fp)

        visualization = Visualization.load(fp)

        self.assertEqual(visualization.type, qiime.core.type.Visualization)
        self.assertEqual(visualization.provenance, self.provenance)
        self.assertEqual(visualization.uuid, saved_visualization.uuid)

    def test_load_and_save(self):
        fp1 = os.path.join(self.test_dir.name, 'visualization1.qzv')
        fp2 = os.path.join(self.test_dir.name, 'visualization2.qzv')
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)
        visualization.save(fp1)

        visualization = Visualization.load(fp1)
        # Overwriting its source file works.
        visualization.save(fp1)
        # Saving to a new file works.
        visualization.save(fp2)

        with zipfile.ZipFile(fp1, mode='r') as zf:
            fps = set(zf.namelist())
            expected = {
                'visualization1/VERSION',
                'visualization1/metadata.yaml',
                'visualization1/README.md',
                'visualization1/data/index.html',
                'visualization1/data/css/style.css'
            }
            self.assertEqual(fps, expected)

        with zipfile.ZipFile(fp2, mode='r') as zf:
            fps = set(zf.namelist())
            expected = {
                'visualization2/VERSION',
                'visualization2/metadata.yaml',
                'visualization2/README.md',
                'visualization2/data/index.html',
                'visualization2/data/css/style.css'
            }
            self.assertEqual(fps, expected)

    def test_roundtrip(self):
        fp1 = os.path.join(self.test_dir.name, 'visualization1.qzv')
        fp2 = os.path.join(self.test_dir.name, 'visualization2.qzv')
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)
        visualization.save(fp1)

        visualization1 = Visualization.load(fp1)
        visualization1.save(fp2)
        visualization2 = Visualization.load(fp2)

        self.assertEqual(visualization1.type, visualization2.type)
        self.assertEqual(visualization1.provenance, visualization2.provenance)
        self.assertEqual(visualization1.uuid, visualization2.uuid)

    def test_load_from_externally_created_zipfile(self):
        # If a user unzips a .qzv to inspect contents and rezips using a
        # different ZIP library/implementation than the one provided by Python,
        # loading, saving, etc. should still work as expected. The Python ZIP
        # implementation doesn't store directories as entries when writing, but
        # the `zip` Unix and OS X command line utilities include both
        # directories and filepaths as entries. When reading these files with
        # Python's ZIP implementation, the directory entries are visible, so
        # their presence needs to be accounted for when extracting.
        #
        # The following visualization was created with:
        #
        #     visualization = Visualization._from_data_dir(self.data_dir,
        #                                                  self.provenance)
        #     visualization.save('externally_created_zipfile.qzv')
        #
        # Unzip and rezip using command line utility:
        #
        #     unzip externally_created_zipfile.qzv
        #     rm externally_created_zipfile.qzv
        #     zip -r externally_created_zipfile.qzv externally_created_zipfile
        #
        fp = pkg_resources.resource_filename(
            'qiime.sdk.tests', 'data/externally_created_zipfile.qzv')

        with zipfile.ZipFile(fp, mode='r') as zf:
            fps = set(zf.namelist())
            expected = {
                # These are extra directory entries included by `zip` command
                # line utility.
                'externally_created_zipfile/',
                'externally_created_zipfile/data/',
                'externally_created_zipfile/data/css/',

                'externally_created_zipfile/VERSION',
                'externally_created_zipfile/metadata.yaml',
                'externally_created_zipfile/README.md',
                'externally_created_zipfile/data/index.html',
                'externally_created_zipfile/data/css/style.css'
            }
            self.assertEqual(fps, expected)

        visualization = Visualization.load(fp)

        self.assertEqual(visualization.type, qiime.core.type.Visualization)
        self.assertEqual(visualization.provenance, self.provenance)
        self.assertIsInstance(visualization.uuid, uuid.UUID)

        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization.save(fp)

        with zipfile.ZipFile(fp, mode='r') as zf:
            fps = set(zf.namelist())
            expected = {
                # Directory entries should not be present.
                'visualization/VERSION',
                'visualization/metadata.yaml',
                'visualization/README.md',
                'visualization/data/index.html',
                'visualization/data/css/style.css'
            }
            self.assertEqual(fps, expected)

    def test_load_with_archive_filepath_modified(self):
        # Save a visualization for use in the following test case.
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        Visualization._from_data_dir(self.data_dir, None).save(fp)

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

        Visualization._from_data_dir(new_data_dir, None).save(fp)
        visualization2 = Visualization.load(fp)

        self.assertEqual(visualization1.get_index_paths(),
                         {'html': 'data/index.html'})
        self.assertEqual(visualization2.get_index_paths(),
                         {'html': 'data/index.html', 'tsv': 'data/index.tsv'})

    def test_extract(self):
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)
        visualization.save(fp)

        output_dir = os.path.join(self.test_dir.name, 'viz-extract-test')
        result_dir = Visualization.extract(fp, output_dir=output_dir)
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

    def test_get_index_paths_single_load(self):
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)
        visualization.save(fp)
        visualization = Visualization.load(fp)

        actual = visualization.get_index_paths()
        expected = {'html': 'data/index.html'}
        self.assertEqual(actual, expected)

    def test_get_index_paths_single_from_data_dir(self):
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)

        actual = visualization.get_index_paths()
        expected = {'html': 'data/index.html'}
        self.assertEqual(actual, expected)

    def test_get_index_paths_multiple_load(self):
        data_dir = os.path.join(self.test_dir.name, 'mc-viz-output1')
        os.mkdir(data_dir)
        most_common_viz(data_dir,
                        collections.Counter(range(42)))
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(data_dir,
                                                     self.provenance)
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
        visualization = Visualization._from_data_dir(data_dir,
                                                     self.provenance)

        actual = visualization.get_index_paths()
        expected = {'html': 'data/index.html',
                    'tsv': 'data/index.tsv'}
        self.assertEqual(actual, expected)

    def test_get_index_paths_multiple_html_load(self):
        data_dir = os.path.join(self.test_dir.name, 'multi-html-viz1')
        os.mkdir(data_dir)
        multi_html_viz(data_dir, [1, 42])

        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization = Visualization._from_data_dir(data_dir,
                                                     self.provenance)
        visualization.save(fp)
        visualization = Visualization.load(fp)

        with self.assertRaises(ValueError):
            visualization.get_index_paths()

    def test_get_index_paths_multiple_html_from_data_dir(self):
        data_dir = os.path.join(self.test_dir.name, 'multi-html-viz2')
        os.mkdir(data_dir)
        multi_html_viz(data_dir, [1, 42])

        visualization = Visualization._from_data_dir(data_dir,
                                                     self.provenance)

        with self.assertRaises(ValueError):
            visualization.get_index_paths()

    def test_get_index_paths_relative_false(self):
        data_dir = os.path.join(self.test_dir.name, 'mc-viz-output2')
        os.mkdir(data_dir)
        most_common_viz(data_dir, collections.Counter(range(42)))
        visualization = Visualization._from_data_dir(data_dir,
                                                     self.provenance)

        def get_abs_path(rel):
            return os.path.join(visualization._archiver._temp_dir, rel)
        actual = visualization.get_index_paths(relative=False)
        expected = {'html': get_abs_path('data/index.html'),
                    'tsv': get_abs_path('data/index.tsv')}
        self.assertEqual(actual, expected)

    def test_peek(self):
        visualization = Visualization._from_data_dir(self.data_dir, None)
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization.save(fp)

        metadata = Visualization.peek(fp)

        self.assertIsInstance(metadata, ResultMetadata)
        self.assertEqual(metadata.type, qiime.core.type.Visualization)
        self.assertIsNone(metadata.provenance)
        self.assertEqual(metadata.uuid, visualization.uuid)

    def test_peek_with_provenance(self):
        visualization = Visualization._from_data_dir(self.data_dir,
                                                     self.provenance)
        fp = os.path.join(self.test_dir.name, 'visualization.qzv')
        visualization.save(fp)

        metadata = Visualization.peek(fp)

        self.assertIsInstance(metadata, ResultMetadata)
        self.assertEqual(metadata.type, qiime.core.type.Visualization)
        self.assertEqual(metadata.provenance, self.provenance)
        self.assertEqual(metadata.uuid, visualization.uuid)


if __name__ == '__main__':
    unittest.main()
