# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import qiime2
import qiime2.sdk


class TestUtil(unittest.TestCase):
    def test_artifact_actions(self):
        obs = qiime2.sdk.util.artifact_actions(None)
        self.assertEqual(obs, [])

        obs = qiime2.sdk.util.artifact_actions('ShouldBeEmpty')
        self.assertEqual(obs, [])

        # For simplicity, we are gonna test the names of the plugin and
        # the actions
        obs = [(x.name, [yy.name for yy in y])
               for x, y in qiime2.sdk.util.artifact_actions('SingleInt')]
        exp = [('dummy-plugin', [
            'Do stuff normally, but override this one step sometimes'])]
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    unittest.main()
