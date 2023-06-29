# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2.sdk.util import view_collection

from qiime2.core.testing.util import get_dummy_plugin


class TestParallel(unittest.TestCase):
    plugin = get_dummy_plugin()

    def test_de_facto_list_arg(self):
        pipeline = self.plugin.pipelines['de_facto_list_pipeline']

        exp = {'0': 0, '1': 1, '2': 2}

        ret = pipeline()
        obs = view_collection(ret.output, int)

        self.assertEqual(obs, exp)

    def test_de_facto_list_arg_parallel(self):
        pipeline = self.plugin.pipelines['de_facto_list_pipeline']

        exp = {'0': 0, '1': 1, '2': 2}

        ret = pipeline.parallel()._result()
        obs = view_collection(ret.output, int)

        self.assertEqual(obs, exp)

    def test_de_facto_list_kwarg(self):
        pipeline = self.plugin.pipelines['de_facto_list_pipeline']

        exp = {'0': 0, '1': 1, '2': 2}

        ret = pipeline(kwarg=True)
        obs = view_collection(ret.output, int)

        self.assertEqual(obs, exp)

    def test_de_facto_list_kwarg_parallel(self):
        pipeline = self.plugin.pipelines['de_facto_list_pipeline']

        exp = {'0': 0, '1': 1, '2': 2}

        ret = pipeline.parallel(kwarg=True)._result()
        obs = view_collection(ret.output, int)

        self.assertEqual(obs, exp)

    def test_de_facto_dict_arg(self):
        pipeline = self.plugin.pipelines['de_facto_dict_pipeline']

        exp = {'1': 0, '2': 1, '3': 2}

        ret = pipeline()
        obs = view_collection(ret.output, int)

        self.assertEqual(obs, exp)

    def test_de_facto_dict_arg_parallel(self):
        pipeline = self.plugin.pipelines['de_facto_dict_pipeline']

        exp = {'1': 0, '2': 1, '3': 2}

        ret = pipeline.parallel()._result()
        obs = view_collection(ret.output, int)

        self.assertEqual(obs, exp)

    def test_de_facto_dict_kwargl(self):
        pipeline = self.plugin.pipelines['de_facto_dict_pipeline']

        exp = {'1': 0, '2': 1, '3': 2}

        ret = pipeline(kwarg=True)
        obs = view_collection(ret.output, int)

        self.assertEqual(obs, exp)

    def test_de_facto_dict_kwarg_parallel(self):
        pipeline = self.plugin.pipelines['de_facto_dict_pipeline']

        exp = {'1': 0, '2': 1, '3': 2}

        ret = pipeline.parallel(kwarg=True)._result()
        obs = view_collection(ret.output, int)

        self.assertEqual(obs, exp)

    def test_fail_de_facto_list_arg_mixed(self):
        pipeline = self.plugin.pipelines['de_facto_list_pipeline']

        with self.assertRaisesRegex(
                ValueError, 'Collection has mixed proxies and artifacts.*'):
            pipeline.parallel(non_proxies=True)._result()

    def test_fail_de_facto_list_kwarg_mixed(self):
        pipeline = self.plugin.pipelines['de_facto_list_pipeline']

        with self.assertRaisesRegex(
                ValueError, 'Collection has mixed proxies and artifacts.*'):
            pipeline.parallel(kwarg=True, non_proxies=True)._result()

    def test_fail_de_facto_dict_arg_mixed(self):
        pipeline = self.plugin.pipelines['de_facto_dict_pipeline']

        with self.assertRaisesRegex(
                ValueError, 'Collection has mixed proxies and artifacts.*'):
            pipeline.parallel(non_proxies=True)._result()

    def test_fail_de_facto_dict_kwarg_mixed(self):
        pipeline = self.plugin.pipelines['de_facto_dict_pipeline']

        with self.assertRaisesRegex(
                ValueError, 'Collection has mixed proxies and artifacts.*'):
            pipeline.parallel(kwarg=True, non_proxies=True)._result()
