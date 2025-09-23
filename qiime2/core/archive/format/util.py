# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pathlib

from contextlib import contextmanager

from qiime2.core.archive import Archiver
from qiime2.core.util import checksum_directory, to_checksum_format


@contextmanager
def artifact_version(version):
    version = str(version)
    if version not in Archiver._FORMAT_REGISTRY:
        raise ValueError("Version %s not supported" % version)
    original_version = Archiver.CURRENT_FORMAT_VERSION
    try:
        Archiver.CURRENT_FORMAT_VERSION = version
        yield
    finally:
        Archiver.CURRENT_FORMAT_VERSION = original_version


def write_checksums(directory, checksum_file, checksum_type):
    checksums = checksum_directory(directory, checksum_type)

    with (pathlib.Path(directory) / checksum_file).open('w') as fh:
        for item in checksums.items():
            # always ignore annotations dir when writing checksums
            # so we dont self-invalidate when adding/removing annotations
            if not item[0].startswith('annotations'):
                fh.write(to_checksum_format(*item))
                fh.write('\n')
