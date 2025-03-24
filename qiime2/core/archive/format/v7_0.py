# ----------------------------------------------------------------------------
# Copyright (c) 2016-2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pathlib

import qiime2.core.archive.format.v1 as v1
import qiime2.core.archive.format.v6 as v6


class ArchiveFormat(v6.ArchiveFormat):
    CONDA_ENV_FILE = 'conda-env.yaml'
    # TODO: NEW annotations dir
    # This will live under prov and can contain license, notes, etc
    # Contains its own self-signed sha256checksum
    # Also references the machine generated md5checksum from top level
    # For any notes files (something like note-1.txt, note-2.txt, etc)
    # Author is required, and this can either be added manually to each note
    # or can be added to the QIIME 2 config with a 'pull default author' flag
    # that can be enabled

    @classmethod
    def write(cls, archive_record, type, format,
              data_initializer, provenance_capture):
        # pulling from the most recent write version that doesn't include
        # checksums - that way we ensure those are only calculated once
        # after all requisite files are present.
        v1.ArchiveFormat.write(archive_record, type, format,
                               data_initializer, provenance_capture)

        # now we add extras within prov specific to v7
        # conda-env.yaml
        conda_fp = \
            archive_record.root / cls.PROVENANCE_DIR / cls.CONDA_ENV_FILE
        conda_prefix = os.environ.get('CONDA_PREFIX')

        if conda_prefix:
            conda_meta_dir = pathlib.Path(conda_prefix) / 'conda-meta'

            if conda_meta_dir.exists():
                meta_files = \
                    [file.stem for file in conda_meta_dir.iterdir()
                     if file.is_file()]

                with conda_fp.open(mode='w') as fh:
                    fh.write('dependencies:\n')
                    fh.writelines(f'- {filename}\n'
                                  for filename in sorted(meta_files))

        else:
            with conda_fp.open(mode='w') as fh:
                fh.write('error: no conda environment detected.\n')

        # md_fp = archive_record.root / cls.METADATA_FILE
        raise ValueError(cls.load_metadata(archive_record))
        # make sure checksums are written last
        cls.write_checksums(archive_record)

    # TODO: figure out how to separate checksum type by self-signed vs.
    # machine generated to ensure that we use sha256 for all self-signed
    # checksums, while all machine generated checksums can remain md5
    @classmethod
    def write_checksums(cls, archive_record):
        super().write_checksums(archive_record)
        # now we write sha256 instead of md5
