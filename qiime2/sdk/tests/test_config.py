# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest
import pkg_resources

import parsl
from parsl.executors.threads import ThreadPoolExecutor
from parsl.errors import NoDataFlowKernelError

from qiime2 import Artifact, Cache
from qiime2.core.util import load_action_yaml
from qiime2.core.testing.type import SingleInt
from qiime2.core.testing.util import get_dummy_plugin
from qiime2.sdk.parallel_config import (PARALLEL_CONFIG, _TEST_EXECUTOR_,
                                        _MASK_CONDA_ENV_, ParallelConfig,
                                        load_config_from_file)


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
        # Create config
        self.test_default = parsl.Config(
            executors=[
                ThreadPoolExecutor(
                    max_threads=1,
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

        # Create temp test dir and cache in dir
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        self.cache = Cache(os.path.join(self.test_dir.name, 'new_cache'))

        # Create artifacts here so we have unique inputs in each test
        self.art = [Artifact.import_data(SingleInt, 0),
                    Artifact.import_data(SingleInt, 1)]

        # Get paths to config files
        self.config_fp = self.get_data_path('test_config.toml')
        self.mapping_config_fp = self.get_data_path('mapping_config.toml')
        self.mapping_only_config_fp = \
            self.get_data_path('mapping_only_config.toml')
        self.complex_config_fp = \
            self.get_data_path('complex_config.toml')

    def tearDown(self):
        self.test_dir.cleanup()

        # Ensure default state post test
        PARALLEL_CONFIG.parallel_config = None
        PARALLEL_CONFIG.action_executor_mapping = {}

    def get_data_path(self, filename):
        return pkg_resources.resource_filename('qiime2.sdk.tests',
                                               'data/%s' % filename)

    def test_default_config(self):
        with ParallelConfig():
            self.assertIsInstance(
                PARALLEL_CONFIG.parallel_config, parsl.Config)
            self.assertEqual(PARALLEL_CONFIG.action_executor_mapping, {})

    def test_mapping_from_config(self):
        config, mapping = load_config_from_file(self.mapping_config_fp)

        with self.cache:
            with ParallelConfig(config, mapping):
                future = self.pipeline.parallel(self.art, self.art)
                list_return, dict_return = future._result()

        list_execution_contexts = self._load_alias_execution_contexts(
            list_return)
        dict_execution_contexts = self._load_alias_execution_contexts(
            dict_return)

        self.assertEqual(list_execution_contexts, self.test_expected)
        self.assertEqual(dict_execution_contexts, self.tpool_expected)

    def test_mapping_only_config(self):
        _, mapping = load_config_from_file(self.mapping_only_config_fp)

        with self.cache:
            with ParallelConfig(action_executor_mapping=mapping):
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
            with ParallelConfig(action_executor_mapping=mapping):
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
                    with ParallelConfig(self.test_default):
                        pass

    def test_parallel_non_pipeline(self):
        with self.assertRaisesRegex(
                ValueError, 'Only pipelines may be run in parallel'):
            self.method.parallel(self.art)

    def test_no_vendored_fp(self):
        with _MASK_CONDA_ENV_():
            with ParallelConfig():
                with self.cache:
                    future = self.pipeline.parallel(self.art, self.art)
                    list_return, dict_return = future._result()

            list_execution_contexts = self._load_alias_execution_contexts(
                list_return)
            dict_execution_contexts = self._load_alias_execution_contexts(
                dict_return)

            self.assertEqual(list_execution_contexts, self.tpool_expected)
            self.assertEqual(dict_execution_contexts, self.tpool_expected)

    def test_load_complex_config(self):
        """ Test that all parsl modules we currently map are correct
        """
        config, mapping = load_config_from_file(self.complex_config_fp)
        # Just assert that we were able to parse the file and get a config out
        with ParallelConfig(config, mapping):
            self.assertIsInstance(
                PARALLEL_CONFIG.parallel_config, parsl.Config)
            self.assertEqual(PARALLEL_CONFIG.action_executor_mapping, {})

    def test_no_config(self):
        with self.assertRaisesRegex(NoDataFlowKernelError,
                                    'Must first load config'):
            self.pipeline.parallel(self.art, self.art)

    def test_config_unset(self):
        with ParallelConfig():
            self.pipeline.parallel(self.art, self.art)

        with self.assertRaisesRegex(NoDataFlowKernelError,
                                    'Must first load config'):
            self.pipeline.parallel(self.art, self.art)

    def _load_alias_execution_contexts(self, collection):
        execution_contexts = []

        for result in collection.values():
            alias_uuid = load_action_yaml(
                result._archiver.path)['action']['alias-of']
            execution_contexts.append(load_action_yaml(
                self.cache.data / alias_uuid)
                ['execution']['execution_context'])

        return execution_contexts
