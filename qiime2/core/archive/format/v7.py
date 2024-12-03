# ----------------------------------------------------------------------------
# Copyright (c) 2016-2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.core.archive.format.v6 as v6


class ArchiveFormat(v6.ArchiveFormat):
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

    # TODO: NEW conda env file
    # This will dump the contents of conda env export into a file
    # that will live under prov/action (something like conda-env.yaml)

    # TODO: UPDATE action.yaml with CPU flags
    # Use psutil to pull these and add as a new section under action.yaml
    pass
