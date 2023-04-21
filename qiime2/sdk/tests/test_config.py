# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import pkg_resources

from parsl import Config

from qiime2.sdk.parsl_config import PARSL_CONFIG, setup_parsl


class TestConfig(unittest.TestCase):
    def setUp(self):
        # Esnure default state prior to test
        PARSL_CONFIG.parsl_config = None
        PARSL_CONFIG.action_executor_mapping = {}

    def get_data_path(self, filename):
        return pkg_resources.resource_filename('qiime2.sdk.tests',
                                               'data/%s' % filename)

    def test_default_config(self):
        config_fp = self.get_data_path('default_config.toml')

        setup_parsl(config_fp)

        # Assert modified state
        self.assertIsInstance(PARSL_CONFIG.parsl_config, Config)
        self.assertEqual(PARSL_CONFIG.action_executor_mapping, {})

    def test_mapping_config(self):
        config_fp = self.get_data_path('mapping_config.toml')

        setup_parsl(config_fp)

        # Assert modified state
        self.assertIsInstance(PARSL_CONFIG.parsl_config, Config)
        self.assertEqual(
            PARSL_CONFIG.action_executor_mapping, {'rarefy': 'tpool'})
