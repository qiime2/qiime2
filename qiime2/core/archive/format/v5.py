# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.core.archive.format.v4 as v4
from qiime2.core.archive.format.util import write_checksums


class ArchiveFormat(v4.ArchiveFormat):
    CHECKSUM_FILE = 'checksums.md5'
    CHECKSUM_TYPE = CHECKSUM_FILE.split('.')[1]

    @classmethod
    def write(cls, archive_record, type, format,
              data_initializer, provenance_capture):
        super().write(archive_record, type, format,
                      data_initializer, provenance_capture)

        # make sure checksums are written last
        cls.write_checksums(archive_record)

    @classmethod
    def write_checksums(cls, archive_record):
        write_checksums(
            directory=str(archive_record.root),
            checksum_file=cls.CHECKSUM_FILE,
            checksum_type=cls.CHECKSUM_TYPE
        )
