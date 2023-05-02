from contextlib import contextmanager
from typing import Generator
import pathlib
import tempfile
import unittest
import zipfile

from ..util import UUID


class CustomAssertions(unittest.TestCase):
    def assertREAppearsOnlyOnce(self, text, only_once, msg=None):
        appears_once_re = \
            (f'(?s)^(?:(?!{only_once}).)*{only_once}(?!.*{only_once}).*$')
        self.assertRegex(text, appears_once_re, msg)


def is_root_provnode_data(fp):
    """
    a filter predicate which returns metadata, action, citation,
    and VERSION fps with which we can construct a ProvNode
    """
    # Handle provenance files...
    if 'provenance' in fp and 'artifacts' not in fp \
        and ('action.yaml' in fp or
             'citations.bib' in fp
             ):
        return True

    # then handle files available at root, which require a cast
    if pathlib.Path(fp).parts[1] in ('VERSION',
                                     'metadata.yaml',
                                     'checksums.md5'):
        return True


@contextmanager
def generate_archive_with_file_removed(qzv_fp: str, root_uuid: UUID,
                                       file_to_drop: pathlib.Path) -> \
                                           Generator[pathlib.Path, None, None]:
    """
    Deleting files from zip archives is hard, so this makes a temporary
    copy of qzf_fp with fp_to_drop removed and returns a handle to this archive

    file_to_drop should represent the relative path to the file within the
    zip archive, excluding the root directory (named for the root UUID).

    e.g. `/d9e080bb-e245-4ab0-a2cf-0a89b63b8050/metadata.yaml` should be passed
    in as `metadata.yaml`

    adapted from https://stackoverflow.com/a/513889/9872253
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_arc = pathlib.Path(tmpdir) / 'mangled.qzv'
        fp_pfx = pathlib.Path(root_uuid)
        zin = zipfile.ZipFile(qzv_fp, 'r')
        zout = zipfile.ZipFile(str(tmp_arc), 'w')
        for item in zin.infolist():
            buffer = zin.read(item.filename)
            drop_filename = str(fp_pfx / file_to_drop)
            if (item.filename != drop_filename):
                zout.writestr(item, buffer)
        zout.close()
        zin.close()
        yield tmp_arc


class ReallyEqualMixin(object):
    """
    Mixin for testing implementations of __eq__/__ne__.

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
