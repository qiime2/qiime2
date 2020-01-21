# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.core.archive.format.v0 as v0


class ArchiveFormat(v0.ArchiveFormat):
    PROVENANCE_DIR = 'provenance'

    @classmethod
    def write(cls, archive_record, type, format, data_initializer,
              provenance_capture):
        super().write(archive_record, type, format, data_initializer,
                      provenance_capture)
        root = archive_record.root

        prov_dir = root / cls.PROVENANCE_DIR
        prov_dir.mkdir()

        provenance_capture.finalize(
            prov_dir, [root / cls.METADATA_FILE, archive_record.version_fp])

    def __init__(self, archive_record):
        super().__init__(archive_record)

        self.provenance_dir = archive_record.root / self.PROVENANCE_DIR
