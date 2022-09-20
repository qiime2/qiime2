# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import shutil
import unittest

from qiime2.plugin import model

# Define dummy plugin formats to test with

class BothRequiredDir(model.DirectoryFormat):
    req1 = model.File(
            'test_text1.txt',
            format=model.TextFileFormat,
            optional=False)

    req2 = model.File(
            'test_text2.txt',
            format=model.TextFileFormat,
            optional=False
            )


class TestDirectoryFormat(unittest.TestCase):
    package = 'qiime2.plugin.model.tests'

    def setUp(self):
        pass

    def testPassingWithAllSomeOptional(self):
        pass

    def testPassingSomeOptional(self):
        pass

    def testPassingAllRequired(self):
        format_object = BothRequiredDir()

        file_one = self.get_data_path('test_text1.txt')
        file_two = self.get_data_path('test_text2.txt')

    def testFailsOnMissingRequired(self):
        pass

