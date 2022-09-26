# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import unittest

import pkg_resources
from qiime2.plugin import model

from qiime2.core.testing.format import IntSequenceFormat
from qiime2.core.exceptions import ValidationError


# Define dummy plugin formats to test with

class AllRequiredDirFmt(model.DirectoryFormat):
    file1 = model.File(r'test_text1.txt', format=IntSequenceFormat,
                       optional=False)
    file2 = model.File(r'test_text2.txt', format=IntSequenceFormat,
                       optional=False)
    file3 = model.File(r'test_text3.txt', format=IntSequenceFormat,
                       optional=False)


class AllRequiredDefaultDirFmt(model.DirectoryFormat):
    file1 = model.File(r'test_text1.txt', format=IntSequenceFormat)
    file2 = model.File(r'test_text2.txt', format=IntSequenceFormat)
    file3 = model.File(r'test_text3.txt', format=IntSequenceFormat)


class OptionalDirFmt(model.DirectoryFormat):
    file1 = model.File(r'test_text1.txt', format=IntSequenceFormat,
                       optional=False)
    file2 = model.File(r'test_text2.txt', format=IntSequenceFormat,
                       optional=False)
    file3 = model.File(r'test_text3.txt', format=IntSequenceFormat,
                       optional=True)


class TestDirectoryFormat(unittest.TestCase):
    package = 'qiime2.plugin.model.tests'

    def get_data_path(self, filename):
        """Convenience method for getting a data asset while testing.

        Test data stored in the ``data/`` dir local to the running test
        can be accessed via this method.

        Parameters
        ----------
        filename : str
            The name of the file to look up.

        Returns
        -------
        filepath : str
            The materialized filepath to the requested test data.

        """

        return pkg_resources.resource_filename(self.package,
                                               'data/%s' % filename)

    def test_fails_missing_required(self):
        files_dir_fp = self.get_data_path('test_text_files/')

        with self.assertRaisesRegex(
            ValidationError, "Missing one or more files for"
                             " AllRequiredDirFmt"):

            format_object = AllRequiredDirFmt(
                                files_dir_fp,
                                mode='r',
                                )

            format_object.validate()

    def test_fails_missing_with_optional_default(self):
        files_dir_fp = self.get_data_path('test_text_files/')

        with self.assertRaisesRegex(ValidationError,
                                    "Missing one or more files for "
                                    "AllRequiredDefaultDirFmt"):
            format_object = AllRequiredDefaultDirFmt(
                                files_dir_fp,
                                mode='r',
                                )
            format_object.validate()

    def test_passes_with_missing_optional(self):
        files_dir_fp = self.get_data_path('test_text_files/')

        format_object = OptionalDirFmt(
                            files_dir_fp,
                            mode='r',
                            )

        format_object.validate()

    def test_fails_on_unknown_file(self):
        files_dir_fp = self.get_data_path('test_text_files_extra/')
        with self.assertRaisesRegex(ValidationError,
                                    ".*Unrecognized file.*"):

            format_object = AllRequiredDirFmt(
                                files_dir_fp,
                                mode='r',
                                )
            format_object.validate()
