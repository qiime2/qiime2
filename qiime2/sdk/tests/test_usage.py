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

    def assert_scope_records_match(self, expected, usage):
        for exp, obs in itertools.zip_longest(expected, usage.scope):
            self.assertEqual(exp.name, obs.name)
            if obs.is_bound:
                self.assertEqual(exp.data_type, str(obs.factory().type))
            else:
                self.assertEqual(exp.data_type, obs.factory)

    def test_no_op_usage_simple(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.NoOpUsage()
        # TODO: looking up actions should be better, probably just a dict
        action.examples[0](use)

        # NoOpUsage assembles records that represent all of the example's
        # provided data (nothing is actually computed)
        expected = [
            X(name='ints_a', data_type='IntSequence1'),
            X(name='ints_b', data_type='IntSequence1'),
            X(name='ints_c', data_type='IntSequence2'),
            X(name='ints_d', data_type=None),
        ]

        self.assert_scope_records_match(expected, use)

    def test_no_op_usage_complex(self):
        action = self.plugin.actions['concatenate_ints']
        use = usage.NoOpUsage()
        action.examples[1](use)

        # NoOpUsage assembles records that represent all of the example's
        # provided data (nothing is actually computed)
        expected = [
            X(name='ints_a', data_type='IntSequence1'),
            X(name='ints_b', data_type='IntSequence1'),
            X(name='ints_c', data_type='IntSequence2'),
            X(name='ints_d', data_type=None),
            X(name='ints_e', data_type=None),
        ]

        self.assert_scope_records_match(expected, use)


if __name__ == '__main__':
    unittest.main()
