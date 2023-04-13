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

from qiime2.sdk.parsl_config import get_config, process_config


class TestConfig(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def get_data_path(self, filename):
        return pkg_resources.resource_filename('qiime2.sdk.tests',
                                               'data/%s' % filename)

    def test_default_config(self):
        config_fp = self.get_data_path('default_config.toml')
        config_dict = get_config(config_fp)
        config_kwargs = process_config(config_dict)
        print(config_kwargs)

        config = Config(**config_kwargs)
        # Not a very useful assertion right now
        self.assertIsInstance(config, Config)