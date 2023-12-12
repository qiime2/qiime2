# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import pkg_resources

import qiime2
import qiime2.sdk

from qiime2.sdk.util import validate_result_collection_keys


class TestUtil(unittest.TestCase):
    def get_data_path(self, filename):
        return pkg_resources.resource_filename('qiime2.sdk.tests',
                                               'data/%s' % filename)

    def test_artifact_actions(self):
        obs = qiime2.sdk.util.actions_by_input_type(None)
        self.assertEqual(obs, [])

        # For simplicity, we are gonna test the names of the plugin and
        # the actions
        # raise ValueError(qiime2.sdk.util.actions_by_input_type('SingleInt'))
        obs = [(x.name, set([yy.name for yy in y]))
               for x, y in qiime2.sdk.util.actions_by_input_type('SingleInt')]
        exp = [('dummy-plugin', set([
            'To be resumed',
            'Do stuff normally, but override this one step sometimes',
            'Internal fail pipeline',
            'Takes and returns a combination of colletions and non collections'
        ]))]
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
            'Visualize most common integers', 'Inputs with typing.Union',
            'Split sequence of integers in half',
            'Test different ways of failing', 'Optional artifacts method',
            'Do stuff normally, but override this one step sometimes',
            'TypeMatch with list and set params'])]
        self.assertEqual(len(obs), 2)
        self.assertEqual(obs[0][0], exp[0][0])
        self.assertCountEqual(obs[0][1], exp[0][1])

    def test_validate_result_collection_keys_valid(self):

        self.assertEqual(validate_result_collection_keys('a'), None)

        good_keys = ['-', '+', '.', '_', 'a', 'x', 'A', 'X', '0', '9',
                     '90XAxa_.+-']
        self.assertEqual(validate_result_collection_keys(*good_keys), None)

    def test_validate_result_collection_keys_invalid(self):
        with self.assertRaisesRegex(KeyError,
                                    "Invalid.*: @"):
            validate_result_collection_keys('@')

        with self.assertRaisesRegex(KeyError,
                                    "Invalid.*: @, a1@"):
            validate_result_collection_keys('@', 'a1@')

        with self.assertRaisesRegex(KeyError,
                                    "Invalid.*: @, a1@"):
            keys = ['@', 'a1@']
            validate_result_collection_keys(*keys)

        with self.assertRaisesRegex(KeyError,
                                    "Invalid.*: @, a1@"):
            validate_result_collection_keys(
                'good-key', '@', 'a1@')

        bad_keys = ['he llo', ' ', '!', '?']
        for key in bad_keys:
            with self.assertRaisesRegex(KeyError,
                                        f"Invalid.*: {key}"):
                validate_result_collection_keys(key)


if __name__ == '__main__':
    unittest.main()
