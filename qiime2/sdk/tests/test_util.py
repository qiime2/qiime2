# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
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
        obs = qiime2.sdk.util.actions_by_input_type(None)
        self.assertEqual(obs, [])

        # For simplicity, we are gonna test the names of the plugin and
        # the actions
        obs = [(x.name, [yy.name for yy in y])
               for x, y in qiime2.sdk.util.actions_by_input_type('SingleInt')]
        exp = [('dummy-plugin', [
            'Do stuff normally, but override this one step sometimes'])]
        self.assertEqual(obs, exp)

        obs = [(x.name, [yy.name for yy in y])
               for x, y in qiime2.sdk.util.actions_by_input_type(
               'Kennel[Cat]')]
        self.assertEqual(obs, [])

        obs = [(x.name, [yy.name for yy in y])
               for x, y in qiime2.sdk.util.actions_by_input_type(
               'IntSequence1')]
        exp = [('dummy-plugin', [
            'A typical pipeline with the potential to raise an error',
            'Concatenate integers', 'Identity', 'Identity', 'Identity',
            'Do a great many things', 'Identity', 'Identity', 'Identity',
            'Visualize most common integers',
            'Split sequence of integers in half',
            'Test different ways of failing', 'Optional artifacts method',
            'Do stuff normally, but override this one step sometimes'])]
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0][0], exp[0][0])
        self.assertCountEqual(obs[0][1], exp[0][1])


if __name__ == '__main__':
    unittest.main()
