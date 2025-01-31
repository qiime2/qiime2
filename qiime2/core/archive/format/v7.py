# ----------------------------------------------------------------------------
# Copyright (c) 2016-2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pathlib

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

    # TODO: NEW filesizes under execution section of action.yaml
    # this will list the total size of all files under data directory

    # We call init_files first to ensure that all files are written prior to
    # checksums being calculated (relevant for this new conda-env.yaml file)
    @classmethod
    def init_files(cls, archive_record, provenance_capture):
        super().init_files(archive_record, provenance_capture)

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

    # need to add a special write operation to ensure that the contents of the
    # data dir are written prior to the prov dir so that file sizes within
    # data dir can be accurately collected and included in action.yaml
    @classmethod
    def write():
        super().write()

    # TODO: figure out how to separate checksum type by self-signed vs.
    # machine generated to ensure that we use sha256 for all self-signed
    # checksums, while all machine generated checksums can remain md5
    @classmethod
    def write_checksums(cls, archive_record):
        super().write_checksums(archive_record)
        # now we write sha256 instead of md5
