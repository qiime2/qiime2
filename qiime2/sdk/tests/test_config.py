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
from parsl.executors import HighThroughputExecutor

from qiime2 import Artifact, Cache
from qiime2.core.util import load_action_yaml
from qiime2.core.testing.type import SingleInt
from qiime2.core.testing.util import get_dummy_plugin
from qiime2.sdk.parsl_config import PARSL_CONFIG, ParallelConfig, setup_parsl


class TestConfig(unittest.TestCase):
    def setUp(self):
        # Ensure default state prior to test
        PARSL_CONFIG.parsl_config = None
        PARSL_CONFIG.action_executor_mapping = {}

        plugin = get_dummy_plugin()
        self.pipeline = plugin.pipelines['resumable_pipeline']

        # Create temp test dir
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')

        # Create artifact and cache
        self.art = [Artifact.import_data(SingleInt, 0),
                    Artifact.import_data(SingleInt, 1)]
        self.cache = Cache(os.path.join(self.test_dir.name, 'new_cache'))

    def tearDown(self):
        self.test_dir.cleanup()

    def get_data_path(self, filename):
        return pkg_resources.resource_filename('qiime2.sdk.tests',
                                               'data/%s' % filename)

    def test_default_config(self):
        config_fp = self.get_data_path('default_config.toml')

        setup_parsl(config_fp)

        # Assert modified state
        self.assertIsInstance(PARSL_CONFIG.parsl_config, Config)
        self.assertEqual(PARSL_CONFIG.action_executor_mapping, {})

    def test_mapping_from_config(self):
        config_fp = self.get_data_path('mapping_config.toml')

        setup_parsl(config_fp)

        # Assert modified state
        self.assertIsInstance(PARSL_CONFIG.parsl_config, Config)
        self.assertEqual(
            PARSL_CONFIG.action_executor_mapping, {'list_of_ints': 'htex'})

        with self.cache:
            future = self.pipeline.parsl(self.art, self.art)
            list_return, dict_return = future.result()

        list_execution_contexts = self._load_alias_execution_contexts(
            list_return)
        dict_execution_contexts = self._load_alias_execution_contexts(
            dict_return)

        list_expected = [{
            'type': 'parsl',
            'parsl_type':
            "<class 'parsl.executors.high_throughput.executor."
            "HighThroughputExecutor'>"}, {
            'type': 'parsl',
            'parsl_type':
            "<class 'parsl.executors.high_throughput.executor."
            "HighThroughputExecutor'>"}]
        dict_expected = [{
            'type': 'parsl',
            'parsl_type':
            "<class 'parsl.executors.threads.ThreadPoolExecutor'>"}, {
            'type': 'parsl',
            'parsl_type':
            "<class 'parsl.executors.threads.ThreadPoolExecutor'>"}]

        self.assertEqual(list_execution_contexts, list_expected)
        self.assertEqual(dict_execution_contexts, dict_expected)

    def test_mapping_from_dict(self):
        config = Config(
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

        mapping = {'list_of_ints': 'htex'}

        with self.cache:
            with ParallelConfig(config, mapping):
                # Assert modified state
                self.assertIsInstance(PARSL_CONFIG.parsl_config, Config)
                self.assertEqual(
                    PARSL_CONFIG.action_executor_mapping,
                    {'list_of_ints': 'htex'})

                future = self.pipeline.parsl(self.art, self.art)
                list_return, dict_return = future.result()

        list_execution_contexts = self._load_alias_execution_contexts(
            list_return)
        dict_execution_contexts = self._load_alias_execution_contexts(
            dict_return)

        list_expected = [{
            'type': 'parsl',
            'parsl_type':
            "<class 'parsl.executors.high_throughput.executor."
            "HighThroughputExecutor'>"}, {
            'type': 'parsl',
            'parsl_type':
            "<class 'parsl.executors.high_throughput.executor."
            "HighThroughputExecutor'>"}]
        dict_expected = [{
            'type': 'parsl',
            'parsl_type':
            "<class 'parsl.executors.threads.ThreadPoolExecutor'>"}, {
            'type': 'parsl',
            'parsl_type':
            "<class 'parsl.executors.threads.ThreadPoolExecutor'>"}]

        self.assertEqual(list_execution_contexts, list_expected)
        self.assertEqual(dict_execution_contexts, dict_expected)

    def _load_alias_execution_contexts(self, collection):
        execution_contexts = []

        for result in collection.values():
            alias_uuid = load_action_yaml(
                result._archiver.path)['action']['alias-of']
            execution_contexts.append(load_action_yaml(
                self.cache.data / alias_uuid)
                ['execution']['execution_context'])

        return execution_contexts
