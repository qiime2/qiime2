# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from pathlib import Path
import os

import qiime2.core.archive.format.v4 as v4
from qiime2.core.util import to_checksum_format, checksum
from qiime2.sdk.result import ChecksumCache


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
        checksum_cache = ChecksumCache()

        checksums = {}
        for root, dirs, files in os.walk(archive_record.root, topdown=True):
            dirs[:] = sorted([d for d in dirs if not d[0] == '.'])
            for file in sorted(files):
                if file[0] == '.':
                    continue

                path = Path(root) / file
                archive_relative_path = path.relative_to(
                    Path(archive_record.root)
                )

                # always ignore annotations dir when checksumming
                # so we dont self-invalidate when adding/removing annotations
                if archive_relative_path.is_relative_to(Path('annotations')):
                    continue

                # check checksum cache
                checksum_string = None
                if archive_relative_path.is_relative_to(
                    Path("provenance/artifacts")
                ):
                    cached_path = archive_relative_path.relative_to(
                        Path("provenance/artifacts")
                    )
                    checksum_string = checksum_cache.get(cached_path)

                if checksum_string is None:
                    checksum_string = checksum(str(path), cls.CHECKSUM_TYPE)

                checksums[str(archive_relative_path)] = checksum_string

        with (Path(archive_record.root) / cls.CHECKSUM_FILE).open('w') as fh:
            for item in checksums.items():
                fh.write(to_checksum_format(*item))
                fh.write('\n')
