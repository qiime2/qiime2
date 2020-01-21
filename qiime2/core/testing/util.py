# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import os.path
import zipfile

import qiime2.sdk


def get_dummy_plugin():
    plugin_manager = qiime2.sdk.PluginManager()
    if 'dummy-plugin' not in plugin_manager.plugins:
        raise RuntimeError(
            "When running QIIME 2 unit tests, the QIIMETEST environment "
            "variable must be defined so that plugins required by unit tests "
            "are loaded. The value of the QIIMETEST environment variable can "
            "be anything. Example command: QIIMETEST=1 nosetests")
    return plugin_manager.plugins['dummy-plugin']


class ArchiveTestingMixin:
    """Mixin for testing properties of archives created by Archiver."""

    def assertArchiveMembers(self, archive_filepath, root_dir, expected):
        """Assert members are in an archive.

        Parameters
        ----------
        archive_filepath : str or Path
            Filepath to archive whose members will be verified against the
            `expected` members.
        root_dir : str or Path
            Root directory of the archive. Will be prepended to the member
            paths in `expected`. This is useful when the archive's root
            directory is not known ahead of time (e.g. when it is a random
            UUID) and the caller is determining the root directory dynamically.
        expected : set of str
            Set of expected archive members stored as paths relative to
            `root_dir`.

        """
        archive_filepath = str(archive_filepath)
        root_dir = str(root_dir)
        with zipfile.ZipFile(archive_filepath, mode='r') as zf:
            observed = set(zf.namelist())

        # Path separator '/' is hardcoded because paths in the zipfile will
        # always use this separator.
        expected = {root_dir + '/' + member for member in expected}

        self.assertEqual(observed, expected)

    def assertExtractedArchiveMembers(self, extract_dir, root_dir, expected):
        """Assert an archive's members are extracted to a directory.

        Parameters
        ----------
        extract_dir : str or Path
            Path to directory the archive was extracted to.
        root_dir : str or Path
            Root directory of the archive that was extracted to `extract_dir`.
            This is useful when the archive's root directory is not known ahead
            of time (e.g. when it is a random UUID) and the caller is
            determining the root directory dynamically.
        expected : set of str
            Set of expected archive members extracted to `extract_dir`. Stored
            as paths relative to `root_dir`.

        """
        extract_dir = str(extract_dir)
        root_dir = str(root_dir)
        observed = set()
        for root, _, filenames in os.walk(extract_dir):
            for filename in filenames:
                observed.add(os.path.join(root, filename))

        expected = {os.path.join(extract_dir, root_dir, member)
                    for member in expected}

        self.assertEqual(observed, expected)


class ReallyEqualMixin:
    """Mixin for testing implementations of __eq__/__ne__.

    Based on this public domain code (also explains why the mixin is useful):

    https://ludios.org/testing-your-eq-ne-cmp/

    """

    def assertReallyEqual(self, a, b):
        # assertEqual first, because it will have a good message if the
        # assertion fails.
        self.assertEqual(a, b)
        self.assertEqual(b, a)
        self.assertTrue(a == b)
        self.assertTrue(b == a)
        self.assertFalse(a != b)
        self.assertFalse(b != a)

    def assertReallyNotEqual(self, a, b):
        # assertNotEqual first, because it will have a good message if the
        # assertion fails.
        self.assertNotEqual(a, b)
        self.assertNotEqual(b, a)
        self.assertFalse(a == b)
        self.assertFalse(b == a)
        self.assertTrue(a != b)
        self.assertTrue(b != a)
