# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.core.archive.format.v4 as v4
from qiime2.core.util import md5sum_directory, to_checksum_format


class ArchiveFormat(v4.ArchiveFormat):
    CHECKSUM_FILE = 'checksums.md5'
    # Adds `checksums.md5` to root of directory structure

    @classmethod
    def write(cls, archive_record, type, format, data_initializer,
              provenance_capture):
        super().write(archive_record, type, format, data_initializer,
                      provenance_capture)

        checksums = md5sum_directory(str(archive_record.root))
        with (archive_record.root / cls.CHECKSUM_FILE).open('w') as fh:
            for item in checksums.items():
                fh.write(to_checksum_format(*item))
                fh.write('\n')
