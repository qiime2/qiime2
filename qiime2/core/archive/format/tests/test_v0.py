# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
import tempfile
import uuid as _uuid
import pathlib
import io

from qiime2.core.testing.type import IntSequence1
from qiime2.core.testing.format import IntSequenceDirectoryFormat
from qiime2.core.archive.archiver import _ZipArchive, ArchiveRecord
from qiime2.core.archive.format.v0 import ArchiveFormat


class TestArchiveFormat(unittest.TestCase):
    def setUp(self):
        prefix = "qiime2-test-temp-"
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)

    def test_format_metadata(self):
        uuid = _uuid.uuid4()
        with io.StringIO() as fh:
            ArchiveFormat._format_metadata(fh, uuid, IntSequence1,
                                           IntSequenceDirectoryFormat)
            result = fh.getvalue()

        self.assertEqual(result,
                         "uuid: %s\ntype: IntSequence1\nformat: "
                         "IntSequenceDirectoryFormat\n" % uuid)

    def test_format_metadata_none(self):
        uuid = _uuid.uuid4()
        with io.StringIO() as fh:
            ArchiveFormat._format_metadata(fh, uuid, IntSequence1, None)
            result = fh.getvalue()

        self.assertEqual(result,
                         "uuid: %s\ntype: IntSequence1\nformat: null\n" % uuid)

    def test_load_root_dir_metadata_uuid_mismatch(self):
        fp = pathlib.Path(self.temp_dir.name) / 'root-dir-metadata-mismatch'
        fp.mkdir()

        r = _ZipArchive.setup(fp, 'foo', 'bar')
        fake = ArchiveRecord(r.root, r.version_fp,
                             _uuid.uuid4(),  # This will trick the format
                             r.version, r.framework_version)

        ArchiveFormat.write(fake, IntSequence1, IntSequenceDirectoryFormat,
                            lambda x: None, None)
        with self.assertRaisesRegex(
                ValueError, 'root directory must match UUID.*metadata'):
            ArchiveFormat(r)


if __name__ == '__main__':
    unittest.main()
