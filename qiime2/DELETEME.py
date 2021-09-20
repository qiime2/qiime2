# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2 import Artifact


class TestExample(unittest.TestCase):
    def test_example(self):
        import qiime2.sdk

        pm = qiime2.sdk.PluginManager()
        self.dummy = pm.get_plugin(id='dummy_plugin')
        a = Artifact.import_data('Foo', "element 1", view_type=str)
        b = Artifact.import_data('Foo', "element 2", view_type=str)

        print(self.dummy)
        future = self.dummy.actions[
            'constrained_input_visualization'].parsl(a=a, b=b)
        future.result()

    test_example2 = test_example
