# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from contextlib import contextmanager

from qiime2.core.archive import Archiver


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
