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

import pandas as pd

import qiime2
from qiime2.core.cache import Cache
from qiime2.core.testing.type import IntSequence1, SingleInt
from qiime2.core.testing.util import get_dummy_plugin, PipelineError
from qiime2.sdk.result import Artifact
from qiime2.core.util import load_action_yaml


def _load_alias_uuid(result):
    return load_action_yaml(result._archiver.path)['action']['alias-of']


def _load_nested_alias_uuid(result, cache):
    alias_uuid = _load_alias_uuid(result)
    aliased_result = qiime2.sdk.Result.load(
        os.path.join(cache.data, alias_uuid))
    return _load_alias_uuid(aliased_result)


def _load_alias_uuids(collection):
    uuids = []

    for result in collection.values():
        uuids.append(_load_alias_uuid(result))

    return uuids


def _load_nested_alias_uuids(collection, cache):
    alias_uuids = _load_alias_uuids(collection)

    # load_alias_uuids is expecting a dictionary, so just make a dictionary
    # here so it gets what it wants
    alias_results = {}
    for idx, alias_uuid in enumerate(alias_uuids):
        alias_results[idx] = \
            qiime2.sdk.Result.load(os.path.join(cache.data, alias_uuid))

    return _load_alias_uuids(alias_results)


class TestPipelineResumption(unittest.TestCase):
    def setUp(self):
        # Get our pipeline
        self.plugin = get_dummy_plugin()
        self.pipeline = self.plugin.pipelines['resumable_varied_pipeline']
        self.nested_pipeline = \
            self.plugin.pipelines['resumable_nested_varied_pipeline']

        # Create temp test dir
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')

        # Create cache and pool
        self.cache = Cache(os.path.join(self.test_dir.name, 'cache'))
        self.pool = self.cache.create_pool('pool')

        # Create artifacts
        self.ints1 = {'1': Artifact.import_data(SingleInt, 0),
                      '2': Artifact.import_data(SingleInt, 1)}
        self.ints1_2 = {'3': Artifact.import_data(SingleInt, 1),
                        '4': Artifact.import_data(SingleInt, 2)}
        self.ints2 = [Artifact.import_data(IntSequence1, [0, 1, 2]),
                      Artifact.import_data(IntSequence1, [3, 4, 5])]
        self.int1 = Artifact.import_data(SingleInt, 42)
        self.int2 = Artifact.import_data(SingleInt, 43)

        # Create metadata
        df1 = pd.DataFrame({'a': ['1', '2', '3']},
                           index=pd.Index(['0', '1', '2'], name='feature ID'))
        self.md1 = qiime2.Metadata(df1)
        df2 = pd.DataFrame({'b': ['4', '5', '6']},
                           index=pd.Index(['0', '1', '2'], name='feature ID'))
        self.md2 = qiime2.Metadata(df2)

    def tearDown(self):
        """Remove our cache and all that from last test
        """
        self.test_dir.cleanup()

    def test_resumable_pipeline_no_pool(self):
        with self.cache:
            with self.assertRaises(PipelineError) as e:
                self.pipeline(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1,
                    fail=True)

            ints1_uuids, ints2_uuids, int1_uuid, list_uuids, dict_uuids, \
                identity_uuid, viz_uuid = e.exception.uuids

            ints1_ret, ints2_ret, int1_ret, list_ret, dict_ret, identity_ret, \
                viz_ret = self.pipeline(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1)

            complete_ints1_uuids = _load_alias_uuids(ints1_ret)
            complete_ints2_uuids = _load_alias_uuids(ints2_ret)
            complete_int1_uuid = _load_alias_uuid(int1_ret)
            complete_list_uuids = _load_alias_uuids(list_ret)
            complete_dict_uuids = _load_alias_uuids(dict_ret)
            complete_identity_uuid = _load_alias_uuid(identity_ret)
            complete_viz_uuid = _load_alias_uuid(viz_ret)

            # Nothing should have been recycled because we didn't use a pool
            self.assertNotEqual(ints1_uuids, complete_ints1_uuids)
            self.assertNotEqual(ints2_uuids, complete_ints2_uuids)
            self.assertNotEqual(int1_uuid, complete_int1_uuid)
            self.assertNotEqual(list_uuids, complete_list_uuids)
            self.assertNotEqual(dict_uuids, complete_dict_uuids)
            self.assertNotEqual(identity_uuid, complete_identity_uuid)
            self.assertNotEqual(viz_uuid, complete_viz_uuid)

    def test_resumable_pipeline(self):
        with self.pool:
            with self.assertRaises(PipelineError) as e:
                self.pipeline(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1,
                    fail=True)

            ints1_uuids, ints2_uuids, int1_uuid, list_uuids, dict_uuids, \
                identity_uuid, viz_uuid = e.exception.uuids

            ints1_ret, ints2_ret, int1_ret, list_ret, dict_ret, \
                identity_ret, viz_ret = self.pipeline(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1)

            complete_ints1_uuids = _load_alias_uuids(ints1_ret)
            complete_ints2_uuids = _load_alias_uuids(ints2_ret)
            complete_int1_uuid = _load_alias_uuid(int1_ret)
            complete_list_uuids = _load_alias_uuids(list_ret)
            complete_dict_uuids = _load_alias_uuids(dict_ret)
            complete_identity_uuid = _load_alias_uuid(identity_ret)
            complete_viz_uuid = _load_alias_uuid(viz_ret)

            # Assert that the artifacts returned by the completed run of
            # the pipeline are aliases of the artifacts created by the
            # first failed run
            self.assertEqual(ints1_uuids, complete_ints1_uuids)
            self.assertEqual(ints2_uuids, complete_ints2_uuids)
            self.assertEqual(int1_uuid, complete_int1_uuid)
            self.assertEqual(list_uuids, complete_list_uuids)
            self.assertEqual(dict_uuids, complete_dict_uuids)
            self.assertEqual(identity_uuid, complete_identity_uuid)
            self.assertEqual(viz_uuid, complete_viz_uuid)

    def test_resumable_pipeline_parallel(self):
        with self.pool:
            with self.assertRaises(PipelineError) as e:
                future = self.pipeline.parallel(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1,
                    fail=True)
                future._result()

            ints1_uuids, ints2_uuids, int1_uuid, list_uuids, dict_uuids, \
                identity_uuid, viz_uuid = e.exception.uuids

            future = self.pipeline.parallel(
                self.ints1, self.ints2, self.int1, 'Hi', self.md1)
            ints1_ret, ints2_ret, int1_ret, list_ret, dict_ret, \
                identity_ret, viz_ret = future._result()

            complete_ints1_uuids = _load_alias_uuids(ints1_ret)
            complete_ints2_uuids = _load_alias_uuids(ints2_ret)
            complete_int1_uuid = _load_alias_uuid(int1_ret)
            complete_list_uuids = _load_alias_uuids(list_ret)
            complete_dict_uuids = _load_alias_uuids(dict_ret)
            complete_identity_uuid = _load_alias_uuid(identity_ret)
            complete_viz_uuid = _load_alias_uuid(viz_ret)

            # Assert that the artifacts returned by the completed run of
            # the pipeline are aliases of the artifacts created by the
            # first failed run
            self.assertEqual(ints1_uuids, complete_ints1_uuids)
            self.assertEqual(ints2_uuids, complete_ints2_uuids)
            self.assertEqual(int1_uuid, complete_int1_uuid)
            self.assertEqual(list_uuids, complete_list_uuids)
            self.assertEqual(dict_uuids, complete_dict_uuids)
            self.assertEqual(identity_uuid, complete_identity_uuid)
            self.assertEqual(viz_uuid, complete_viz_uuid)

    def test_resumable_pipeline_artifact_varies(self):
        with self.pool:
            with self.assertRaises(PipelineError) as e:
                self.pipeline(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1,
                    fail=True)

            ints1_uuids, ints2_uuids, int1_uuid, list_uuids, dict_uuids, \
                identity_uuid, viz_uuid = e.exception.uuids

            # Pass int2 instead of int1
            ints1_ret, ints2_ret, int1_ret, list_ret, dict_ret, \
                identity_ret, viz_ret = self.pipeline(
                    self.ints1, self.ints2, self.int2, 'Hi', self.md1)

            complete_ints1_uuids = _load_alias_uuids(ints1_ret)
            complete_ints2_uuids = _load_alias_uuids(ints2_ret)
            complete_int1_uuid = _load_alias_uuid(int1_ret)
            complete_list_uuids = _load_alias_uuids(list_ret)
            complete_dict_uuids = _load_alias_uuids(dict_ret)
            complete_identity_uuid = _load_alias_uuid(identity_ret)
            complete_viz_uuid = _load_alias_uuid(viz_ret)

            # Assert that the artifacts returned by the completed pipeline that
            # are implicated by the changed input are not aliases while the
            # others are
            self.assertNotEqual(ints1_uuids, complete_ints1_uuids)
            self.assertNotEqual(ints2_uuids, complete_ints2_uuids)
            self.assertNotEqual(int1_uuid, complete_int1_uuid)
            self.assertNotEqual(list_uuids, complete_list_uuids)
            self.assertEqual(dict_uuids, complete_dict_uuids)
            self.assertEqual(identity_uuid, complete_identity_uuid)
            self.assertEqual(viz_uuid, complete_viz_uuid)

    def test_resumable_pipeline_artifact_varies_parallel(self):
        with self.pool:
            with self.assertRaises(PipelineError) as e:
                future = self.pipeline.parallel(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1,
                    fail=True)
                future._result()

            ints1_uuids, ints2_uuids, int1_uuid, list_uuids, dict_uuids, \
                identity_uuid, viz_uuid = e.exception.uuids

            # Pass int2 instead of int1
            future = self.pipeline.parallel(
                self.ints1, self.ints2, self.int2, 'Hi', self.md1)
            ints1_ret, ints2_ret, int1_ret, list_ret, dict_ret, identity_ret, \
                viz_ret = future._result()

            complete_ints1_uuids = _load_alias_uuids(ints1_ret)
            complete_ints2_uuids = _load_alias_uuids(ints2_ret)
            complete_int1_uuid = _load_alias_uuid(int1_ret)
            complete_list_uuids = _load_alias_uuids(list_ret)
            complete_dict_uuids = _load_alias_uuids(dict_ret)
            complete_identity_uuid = _load_alias_uuid(identity_ret)
            complete_viz_uuid = _load_alias_uuid(viz_ret)

            # Assert that the artifacts returned by the completed pipeline that
            # are implicated by the changed input are not aliases while the
            # others are
            self.assertNotEqual(ints1_uuids, complete_ints1_uuids)
            self.assertNotEqual(ints2_uuids, complete_ints2_uuids)
            self.assertNotEqual(int1_uuid, complete_int1_uuid)
            self.assertNotEqual(list_uuids, complete_list_uuids)
            self.assertEqual(dict_uuids, complete_dict_uuids)
            self.assertEqual(identity_uuid, complete_identity_uuid)
            self.assertEqual(viz_uuid, complete_viz_uuid)

    def test_resumable_pipeline_collection_varies(self):
        with self.pool:
            with self.assertRaises(PipelineError) as e:
                self.pipeline(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1,
                    fail=True)

            ints1_uuids, ints2_uuids, int1_uuid, list_uuids, dict_uuids, \
                identity_uuid, viz_uuid = e.exception.uuids

            # Pass ints1_2 instead of ints1
            ints1_ret, ints2_ret, int1_ret, list_ret, dict_ret, \
                identity_ret, viz_ret = self.pipeline(
                    self.ints1_2, self.ints2, self.int1, 'Hi', self.md1)

            complete_ints1_uuids = _load_alias_uuids(ints1_ret)
            complete_ints2_uuids = _load_alias_uuids(ints2_ret)
            complete_int1_uuid = _load_alias_uuid(int1_ret)
            complete_list_uuids = _load_alias_uuids(list_ret)
            complete_dict_uuids = _load_alias_uuids(dict_ret)
            complete_identity_uuid = _load_alias_uuid(identity_ret)
            complete_viz_uuid = _load_alias_uuid(viz_ret)

            # Assert that the artifacts returned by the completed pipeline that
            # are implicated by the changed input are not aliases while the
            # others are
            self.assertNotEqual(ints1_uuids, complete_ints1_uuids)
            self.assertNotEqual(ints2_uuids, complete_ints2_uuids)
            self.assertNotEqual(int1_uuid, complete_int1_uuid)
            self.assertNotEqual(list_uuids, complete_list_uuids)
            self.assertNotEqual(dict_uuids, complete_dict_uuids)
            self.assertEqual(identity_uuid, complete_identity_uuid)
            self.assertEqual(viz_uuid, complete_viz_uuid)

    def test_resumable_pipeline_collection_varies_parallel(self):
        with self.pool:
            with self.assertRaises(PipelineError) as e:
                future = self.pipeline.parallel(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1,
                    fail=True)
                future._result()

            ints1_uuids, ints2_uuids, int1_uuid, list_uuids, dict_uuids, \
                identity_uuid, viz_uuid = e.exception.uuids

            # Pass ints1_2 instead of ints1
            future = self.pipeline.parallel(
                self.ints1_2, self.ints2, self.int2, 'Hi', self.md1)
            ints1_ret, ints2_ret, int1_ret, list_ret, dict_ret, identity_ret, \
                viz_ret = future._result()

            complete_ints1_uuids = _load_alias_uuids(ints1_ret)
            complete_ints2_uuids = _load_alias_uuids(ints2_ret)
            complete_int1_uuid = _load_alias_uuid(int1_ret)
            complete_list_uuids = _load_alias_uuids(list_ret)
            complete_dict_uuids = _load_alias_uuids(dict_ret)
            complete_identity_uuid = _load_alias_uuid(identity_ret)
            complete_viz_uuid = _load_alias_uuid(viz_ret)

            # Assert that the artifacts returned by the completed pipeline that
            # are implicated by the changed input are not aliases while the
            # others are
            self.assertNotEqual(ints1_uuids, complete_ints1_uuids)
            self.assertNotEqual(ints2_uuids, complete_ints2_uuids)
            self.assertNotEqual(int1_uuid, complete_int1_uuid)
            self.assertNotEqual(list_uuids, complete_list_uuids)
            self.assertNotEqual(dict_uuids, complete_dict_uuids)
            self.assertEqual(identity_uuid, complete_identity_uuid)
            self.assertEqual(viz_uuid, complete_viz_uuid)

    def test_resumable_pipeline_str_varies(self):
        with self.pool:
            with self.assertRaises(PipelineError) as e:
                self.pipeline(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1,
                    fail=True)

            ints1_uuids, ints2_uuids, int1_uuid, list_uuids, dict_uuids, \
                identity_uuid, viz_uuid = e.exception.uuids

            # Pass in Bye instead of Hi
            ints1_ret, ints2_ret, int1_ret, list_ret, dict_ret, \
                identity_ret, viz_ret = self.pipeline(
                    self.ints1, self.ints2, self.int1, 'Bye', self.md1)

            complete_ints1_uuids = _load_alias_uuids(ints1_ret)
            complete_ints2_uuids = _load_alias_uuids(ints2_ret)
            complete_int1_uuid = _load_alias_uuid(int1_ret)
            complete_list_uuids = _load_alias_uuids(list_ret)
            complete_dict_uuids = _load_alias_uuids(dict_ret)
            complete_identity_uuid = _load_alias_uuid(identity_ret)
            complete_viz_uuid = _load_alias_uuid(viz_ret)

            # Assert that the artifacts returned by the completed pipeline that
            # are implicated by the changed input are not aliases while the
            # others are
            self.assertNotEqual(ints1_uuids, complete_ints1_uuids)
            self.assertNotEqual(ints2_uuids, complete_ints2_uuids)
            self.assertNotEqual(int1_uuid, complete_int1_uuid)
            self.assertNotEqual(list_uuids, complete_list_uuids)
            self.assertEqual(dict_uuids, complete_dict_uuids)
            self.assertEqual(identity_uuid, complete_identity_uuid)
            self.assertEqual(viz_uuid, complete_viz_uuid)

    def test_resumable_pipeline_str_varies_parallel(self):
        with self.pool:
            with self.assertRaises(PipelineError) as e:
                future = self.pipeline.parallel(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1,
                    fail=True)
                future._result()

            ints1_uuids, ints2_uuids, int1_uuid, list_uuids, dict_uuids, \
                identity_uuid, viz_uuid = e.exception.uuids

            # Pass in Bye instead of Hi
            future = self.pipeline.parallel(
                self.ints1, self.ints2, self.int1, 'Bye', self.md1)
            ints1_ret, ints2_ret, int1_ret, list_ret, dict_ret, \
                identity_ret, viz_ret = future._result()

            complete_ints1_uuids = _load_alias_uuids(ints1_ret)
            complete_ints2_uuids = _load_alias_uuids(ints2_ret)
            complete_int1_uuid = _load_alias_uuid(int1_ret)
            complete_list_uuids = _load_alias_uuids(list_ret)
            complete_dict_uuids = _load_alias_uuids(dict_ret)
            complete_identity_uuid = _load_alias_uuid(identity_ret)
            complete_viz_uuid = _load_alias_uuid(viz_ret)

            # Assert that the artifacts returned by the completed pipeline that
            # are implicated by the changed input are not aliases while the
            # others are
            self.assertNotEqual(ints1_uuids, complete_ints1_uuids)
            self.assertNotEqual(ints2_uuids, complete_ints2_uuids)
            self.assertNotEqual(int1_uuid, complete_int1_uuid)
            self.assertNotEqual(list_uuids, complete_list_uuids)
            self.assertEqual(dict_uuids, complete_dict_uuids)
            self.assertEqual(identity_uuid, complete_identity_uuid)
            self.assertEqual(viz_uuid, complete_viz_uuid)

    def test_resumable_pipeline_md_varies(self):
        with self.pool:
            with self.assertRaises(PipelineError) as e:
                self.pipeline(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1,
                    fail=True)

            ints1_uuids, ints2_uuids, int1_uuid, list_uuids, dict_uuids, \
                identity_uuid, viz_uuid = e.exception.uuids

            # Pass in md2 instead of md1
            ints1_ret, ints2_ret, int1_ret, list_ret, dict_ret, \
                identity_ret, viz_ret = self.pipeline(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md2)

            complete_ints1_uuids = _load_alias_uuids(ints1_ret)
            complete_ints2_uuids = _load_alias_uuids(ints2_ret)
            complete_int1_uuid = _load_alias_uuid(int1_ret)
            complete_list_uuids = _load_alias_uuids(list_ret)
            complete_dict_uuids = _load_alias_uuids(dict_ret)
            complete_identity_uuid = _load_alias_uuid(identity_ret)
            complete_viz_uuid = _load_alias_uuid(viz_ret)

            # Assert that the artifacts returned by the completed pipeline that
            # are implicated by the changed input are not aliases while the
            # others are
            self.assertEqual(ints1_uuids, complete_ints1_uuids)
            self.assertEqual(ints2_uuids, complete_ints2_uuids)
            self.assertEqual(int1_uuid, complete_int1_uuid)
            self.assertEqual(list_uuids, complete_list_uuids)
            self.assertEqual(dict_uuids, complete_dict_uuids)
            self.assertNotEqual(identity_uuid, complete_identity_uuid)
            self.assertEqual(viz_uuid, complete_viz_uuid)

    def test_resumable_pipeline_md_varies_parallel(self):
        with self.pool:
            with self.assertRaises(PipelineError) as e:
                future = self.pipeline.parallel(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1,
                    fail=True)
                future._result()

            ints1_uuids, ints2_uuids, int1_uuid, list_uuids, dict_uuids, \
                identity_uuid, viz_uuid = e.exception.uuids

            # Pass in md2 instead of md1
            future = self.pipeline.parallel(
                self.ints1, self.ints2, self.int1, 'Hi', self.md2)
            ints1_ret, ints2_ret, int1_ret, list_ret, dict_ret, \
                identity_ret, viz_ret = future._result()

            complete_ints1_uuids = _load_alias_uuids(ints1_ret)
            complete_ints2_uuids = _load_alias_uuids(ints2_ret)
            complete_int1_uuid = _load_alias_uuid(int1_ret)
            complete_list_uuids = _load_alias_uuids(list_ret)
            complete_dict_uuids = _load_alias_uuids(dict_ret)
            complete_identity_uuid = _load_alias_uuid(identity_ret)
            complete_viz_uuid = _load_alias_uuid(viz_ret)

            # Assert that the artifacts returned by the completed pipeline that
            # are implicated by the changed input are not aliases while the
            # others are
            self.assertEqual(ints1_uuids, complete_ints1_uuids)
            self.assertEqual(ints2_uuids, complete_ints2_uuids)
            self.assertEqual(int1_uuid, complete_int1_uuid)
            self.assertEqual(list_uuids, complete_list_uuids)
            self.assertEqual(dict_uuids, complete_dict_uuids)
            self.assertNotEqual(identity_uuid, complete_identity_uuid)
            self.assertEqual(viz_uuid, complete_viz_uuid)

    def test_nested_resumable_pipeline(self):
        with self.pool:
            with self.assertRaises(PipelineError) as e:
                self.nested_pipeline(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1,
                    fail=True)

            ints1_uuids, ints2_uuids, int1_uuid, list_uuids, dict_uuids, \
                identity_uuid, viz_uuid = e.exception.uuids

            # We now run the not nested version. This will be able to reuse the
            # returns from varied_method
            ints1_ret, ints2_ret, int1_ret, list_ret, dict_ret, \
                identity_ret, viz_ret = self.nested_pipeline(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1)

            complete_ints1_uuids = _load_nested_alias_uuids(
                ints1_ret, self.cache)
            complete_ints2_uuids = _load_nested_alias_uuids(
                ints2_ret, self.cache)
            complete_int1_uuid = _load_nested_alias_uuid(int1_ret, self.cache)
            complete_list_uuids = _load_alias_uuids(list_ret)
            complete_dict_uuids = _load_alias_uuids(dict_ret)
            complete_identity_uuid = _load_alias_uuid(identity_ret)
            complete_viz_uuid = _load_alias_uuid(viz_ret)

            # Assert that the artifacts returned by the completed run of the
            # pipeline are aliases of the artifacts created by the first failed
            # run
            self.assertEqual(ints1_uuids, complete_ints1_uuids)
            self.assertEqual(ints2_uuids, complete_ints2_uuids)
            self.assertEqual(int1_uuid, complete_int1_uuid)
            self.assertEqual(list_uuids, complete_list_uuids)
            self.assertEqual(dict_uuids, complete_dict_uuids)
            self.assertEqual(identity_uuid, complete_identity_uuid)
            self.assertEqual(viz_uuid, complete_viz_uuid)

    def test_nested_resumable_pipeline_parallel(self):
        with self.pool:
            with self.assertRaises(PipelineError) as e:
                future = self.nested_pipeline.parallel(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1,
                    fail=True)
                future._result()

            ints1_uuids, ints2_uuids, int1_uuid, list_uuids, dict_uuids, \
                identity_uuid, viz_uuid = e.exception.uuids

            future = self.nested_pipeline.parallel(
                    self.ints1, self.ints2, self.int1, 'Hi', self.md1)
            ints1_ret, ints2_ret, int1_ret, list_ret, dict_ret, \
                identity_ret, viz_ret = future._result()

            complete_ints1_uuids = _load_nested_alias_uuids(
                ints1_ret, self.cache)
            complete_ints2_uuids = _load_nested_alias_uuids(
                ints2_ret, self.cache)
            complete_int1_uuid = _load_nested_alias_uuid(int1_ret, self.cache)
            complete_list_uuids = _load_alias_uuids(list_ret)
            complete_dict_uuids = _load_alias_uuids(dict_ret)
            complete_identity_uuid = _load_alias_uuid(identity_ret)
            complete_viz_uuid = _load_alias_uuid(viz_ret)

            # Assert that the artifacts returned by the completed run of the
            # pipeline are aliases of the artifacts created by the first failed
            # run
            self.assertEqual(ints1_uuids, complete_ints1_uuids)
            self.assertEqual(ints2_uuids, complete_ints2_uuids)
            self.assertEqual(int1_uuid, complete_int1_uuid)
            self.assertEqual(list_uuids, complete_list_uuids)
            self.assertEqual(dict_uuids, complete_dict_uuids)
            self.assertEqual(identity_uuid, complete_identity_uuid)
            self.assertEqual(viz_uuid, complete_viz_uuid)
