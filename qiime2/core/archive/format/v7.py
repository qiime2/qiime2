# ----------------------------------------------------------------------------
# Copyright (c) 2016-2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess

import qiime2.core.archive.format.v6 as v6


class ArchiveFormat(v6.ArchiveFormat):
    CONDA_ENV_FILE = 'conda-env.yaml'
    # TODO: NEW annotations dir
    # This will live under prov and can contain license, notes, etc
    # Contains its own self-signed checksum
    # Also references the machine generated checksum from top level
    # For any notes files (something like note-1.txt, note-2.txt, etc)
    # Author is required, and this can either be added manually to each note
    # or can be added to the QIIME 2 config with a 'pull default author' flag
    # that can be enabled

    # TODO: NEW filesizes file
    # This will list all of the files within the data dir
    # and their respective sizes in bytes
    # This will be its own file in the top level dir
    # (similar to the checksums file)

    # TODO: UPDATE action.yaml with CPU flags
    # Use psutil to pull these and add as a new section under action.yaml

    # We call init_files first to ensure that all files are written prior to
    # checksums being calculated (relevant for this new conda-env.yaml file)
    @classmethod
    def init_files(cls, archive_record, provenance_capture):
        super().init_files(archive_record, provenance_capture)

        conda_fp = \
            archive_record.root / cls.PROVENANCE_DIR / cls.CONDA_ENV_FILE

        try:
            cmd = subprocess.run(["conda", "env", "export"],
                                 capture_output=True,
                                 text=True, check=True)

            lines = cmd.stdout.splitlines()
            filtered_lines = [
                line + '\n' for line in lines
                if not (line.startswith("name:") or line.startswith("prefix:"))
            ]

            with conda_fp.open(mode='w') as fh:
                fh.writelines(filtered_lines)

        except subprocess.CalledProcessError as e:
            print(f"Error exporting conda environment: {e}")

    # Now that all files are written, can write the checksums file
    # for everyone, now using sha256 instead of md5
    @classmethod
    def write_checksums(cls, archive_record):
        super().write_checksums(archive_record)
        # now we write sha256 instead of md5
