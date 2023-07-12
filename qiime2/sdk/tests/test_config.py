# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import psutil
import tempfile
import unittest
import pkg_resources

from parsl import Config
from parsl.providers import LocalProvider
from parsl.executors.threads import ThreadPoolExecutor

from qiime2 import Artifact, Cache
from qiime2.core.util import load_action_yaml
from qiime2.core.testing.type import SingleInt
from qiime2.core.testing.util import get_dummy_plugin
from qiime2.sdk.parallel_config import (PARALLEL_CONFIG, PallelConfig,
                                        setup_parallel)


class TestConfig(unittest.TestCase):
    # Get actions
    plugin = get_dummy_plugin()
    pipeline = plugin.pipelines['resumable_pipeline']
    method = plugin.methods['list_of_ints']

    # Expected provenance based on type of executor used
    tpool_expected = [{
        'type': 'parsl', 'parsl_type': 'ThreadPoolExecutor'}, {
        'type': 'parsl', 'parsl_type': 'ThreadPoolExecutor'}]
    test_expected = [{
        'type': 'parsl', 'parsl_type': '_TEST_EXECUTOR_'}, {
        'type': 'parsl', 'parsl_type': '_TEST_EXECUTOR_'}]

    def setUp(self):
        # Ensure default state prior to test
        PARALLEL_CONFIG.parallel_config = None
        PARALLEL_CONFIG.action_executor_mapping = {}

        self.test_default = parsl.Config(
            executors=[
                ThreadPoolExecutor(
                    max_threads=max(psutil.cpu_count() - 1, 1),
                    label='default'
                ),
                HighThroughputExecutor(
                    label='htex',
                    max_workers=max(psutil.cpu_count() - 1, 1),
                    provider=LocalProvider()
                )
            ],
            # AdHoc Clusters should not be setup with scaling strategy.
            strategy='none',
        )

        self.htex_default = Config(
            executors=[
                ThreadPoolExecutor(
                    max_threads=max(psutil.cpu_count() - 1, 1),
                    label='tpool'
                ),
                _TEST_EXECUTOR_(
                    max_threads=1,
                    label='default'
                )
            ],
            # AdHoc Clusters should not be setup with scaling strategy.
            strategy='none',
        )

        self.tpool_expected = [{
            'type': 'parsl', 'parsl_type': 'ThreadPoolExecutor'}, {
            'type': 'parsl', 'parsl_type': 'ThreadPoolExecutor'}]
        self.htex_expected = [{
            'type': 'parsl', 'parsl_type': 'HighThroughputExecutor'}, {
            'type': 'parsl', 'parsl_type': 'HighThroughputExecutor'}]

    def tearDown(self):
        self.test_dir.cleanup()

    def get_data_path(self, filename):
        return pkg_resources.resource_filename('qiime2.sdk.tests',
                                               'data/%s' % filename)

    def test_default_config(self):
        config_fp = self.get_data_path('default_config.toml')

        setup_parallel(config_fp)

        # Assert modified state
        self.assertIsInstance(PARALLEL_CONFIG.parallel_config, Config)
        self.assertEqual(PARALLEL_CONFIG.action_executor_mapping, {})

    def test_mapping_from_config(self):
        setup_parallel(self.config_fp)

        with self.cache:
            future = self.pipeline.parallel(self.art, self.art)
            list_return, dict_return = future._result()

        list_execution_contexts = self._load_alias_execution_contexts(
            list_return)
        dict_execution_contexts = self._load_alias_execution_contexts(
            dict_return)

        self.assertEqual(list_execution_contexts, self.test_expected)
        self.assertEqual(dict_execution_contexts, self.tpool_expected)

    def test_mapping_from_dict(self):
        mapping = {'list_of_ints': 'test'}

        with self.cache:
            with ParallelConfig(self.config, mapping):
                future = self.pipeline.parallel(self.art, self.art)
                list_return, dict_return = future._result()

        list_execution_contexts = self._load_alias_execution_contexts(
            list_return)
        dict_execution_contexts = self._load_alias_execution_contexts(
            dict_return)

        self.assertEqual(list_execution_contexts, self.test_expected)
        self.assertEqual(dict_execution_contexts, self.tpool_expected)

    def test_parallel_configs(self):
        with self.cache:
            with ParallelConfig():
                future = self.pipeline.parallel(self.art, self.art)
                list_return, dict_return = future._result()

            list_execution_contexts = self._load_alias_execution_contexts(
                list_return)
            dict_execution_contexts = self._load_alias_execution_contexts(
                dict_return)

            self.assertEqual(list_execution_contexts, self.tpool_expected)
            self.assertEqual(dict_execution_contexts, self.tpool_expected)

            with ParallelConfig(self.test_default):
                future = self.pipeline.parallel(self.art, self.art)
                list_return, dict_return = future._result()

            list_execution_contexts = self._load_alias_execution_contexts(
                list_return)
            dict_execution_contexts = self._load_alias_execution_contexts(
                dict_return)

            self.assertEqual(list_execution_contexts, self.test_expected)
            self.assertEqual(dict_execution_contexts, self.test_expected)

            # At this point we should be using the default config again which
            # does not have an executor called tpool
            with ParallelConfig(
                    action_executor_mapping={'list_of_ints': 'tpool'}):
                with self.assertRaisesRegex(KeyError, 'tpool'):
                    future = self.pipeline.parallel(self.art, self.art)
                    list_return, dict_return = future._result()

    def test_nested_configs(self):
        with self.cache:
            with self.assertRaisesRegex(
                    ValueError, 'cannot nest ParallelConfigs'):
                with ParallelConfig():
                    with ParallelConfig():
                        pass

    def test_parallel_non_pipeline(self):
        with self.assertRaisesRegex(
                ValueError, 'Only pipelines may be run in parallel'):
            self.method.parallel(self.art)

    def test_no_vendored_fp(self):
        with _MaskCondaEnv():
            with ParallelConfig():
                with self.cache:
                    future = self.pipeline.parallel(self.art, self.art)
                    list_return, dict_return = future._result()

                    list_execution_contexts = \
                        self._load_alias_execution_contexts(list_return)
                    dict_execution_contexts = \
                        self._load_alias_execution_contexts(dict_return)

                    self.assertEqual(
                        list_execution_contexts, self.htex_expected)
                    self.assertEqual(
                        dict_execution_contexts, self.htex_expected)

                future = self.pipeline.parallel(self.art, self.art)
                list_return, dict_return = future._result()

                list_execution_contexts = self._load_alias_execution_contexts(
                    list_return)
                dict_execution_contexts = self._load_alias_execution_contexts(
                    dict_return)

                self.assertEqual(list_execution_contexts, self.tpool_expected)
                self.assertEqual(dict_execution_contexts, self.tpool_expected)

    def test_parallel_non_pipeline(self):
        with self.assertRaisesRegex(
                ValueError, 'Only pipelines may be run in parallel'):
            self.method.parallel(self.art)

    def _load_alias_execution_contexts(self, collection):
        execution_contexts = []

        for result in collection.values():
            alias_uuid = load_action_yaml(
                result._archiver.path)['action']['alias-of']
            execution_contexts.append(load_action_yaml(
                self.cache.data / alias_uuid)
                ['execution']['execution_context'])

        return execution_contexts
