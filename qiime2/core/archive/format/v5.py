# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.core.archive.format.v4 as v4
from qiime2.core.util import md5sum_directory, to_checksum_format


class ArchiveFormat(v4.ArchiveFormat):
    # Adds `checksums.md5` to root of directory structure
    CHECKSUM_FILE = 'checksums.md5'

    @classmethod
    def write(cls, archive_record, type, format,
              data_initializer, provenance_capture):
        super().write(archive_record, type, format,
                      data_initializer, provenance_capture)

        # make sure checksums are written last
        cls.write_checksums(archive_record)

    @classmethod
    def write_checksums(cls, archive_record):
        checksums = md5sum_directory(str(archive_record.root))
        with (archive_record.root / cls.CHECKSUM_FILE).open('w') as fh:
            for item in checksums.items():
                fh.write(to_checksum_format(*item))
                fh.write('\n')
