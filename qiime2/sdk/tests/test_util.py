# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
from qiime2.core.testing.util import get_dummy_plugin
from qiime2 import Artifact
from qiime2.sdk.util import has_metadata


class TestUtil(unittest.TestCase):

    def setUp(self):
        get_dummy_plugin()

    def test_has_metadata_true(self):
        A = Artifact.import_data('Mapping', {'a': '1', 'b': '2'})
        self.assertTrue(has_metadata(A))

    def test_has_metadata_false(self):
        A = Artifact.import_data('IntSequence1', [1, 2, 3, 4])
        self.assertFalse(has_metadata(A))


if __name__ == "__main__":
    unittest.main()
