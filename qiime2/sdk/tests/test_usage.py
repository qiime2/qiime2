# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import itertools
import os
import tempfile
import unittest

from qiime2.sdk import usage
from qiime2.core.testing.util import get_dummy_plugin
from qiime2.plugin.testing import TestPluginBase


X = collections.namedtuple('X', 'name data_type')


class TestTesting(TestPluginBase):
    def setUp(self):
        self.plugin = get_dummy_plugin()

        # TODO standardize temporary directories created by QIIME 2
        # create a temporary data_dir for sample Visualizations
        self.test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        self.data_dir = os.path.join(self.test_dir.name, 'viz-output')
        os.mkdir(self.data_dir)

    def tearDown(self):
        self.test_dir.cleanup()

    def test_no_op_usage(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.NoOpUsage()

        for example in action.examples:
            example(use)

        # NoOpUsage assembles records that represent all of the example's
        # provided data (nothing is actually computed)
        expected = [
            X(name='byod', data_type='IntSequence1'),
            X(name='ints2', data_type='IntSequence1'),
            X(name='this_one_is_important', data_type='IntSequence2'),
            X(name='youre_just_a_copy_of_an_imitation', data_type='None'),
            X(name='well_well_well_what_do_we_have_here', data_type='None'),
        ]

        for exp, obs in itertools.zip_longest(expected, use.scope):
            self.assertEqual(exp.name, obs.name)
            # TODO: fix this
            self.assertEqual(exp.data_type, str(type(obs.factory().type)))


if __name__ == '__main__':
    unittest.main()
