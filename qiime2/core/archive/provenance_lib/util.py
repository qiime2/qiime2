# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import codecs
import pathlib
import re
import warnings

from typing import Tuple
from zipfile import ZipFile


def get_root_uuid(zf: ZipFile) -> str:
    '''
    Returns the root UUID of a QIIME 2 Archive.

    Parameters
    ----------
    zf : ZipFile
        The zipfile object of an archive.

    Returns
    -------
    str
        The uuid of the root artifact in the archive.
    '''
    return pathlib.Path(zf.namelist()[0]).parts[0]


def get_nonroot_uuid(fp: pathlib.Path) -> str:
    '''
    For non-root provenance files, get the Result's uuid from its path.

    Parameters
    ----------
    fp : pathlib.Path
        The path to a file in a non-root artifact inside an archive, relative
        to archive root.

    Returns
    -------
    str
        The uuid of the non-root artifact.
    '''
    if fp.name == 'action.yaml':
        return fp.parts[-3]
    return fp.parts[-2]


_VERSION_MATCHER = (
    r'QIIME 2\n'
    r'archive: [0-9]{1,2}$\n'
    r'framework: '
    r'(?:20[0-9]{2}|2)\.(?:[1-9][0-2]?|0)\.[0-9](?:\.dev[0-9]?)?'
    r'(?:\+[.\w]+)?\Z'
)


def parse_version(zf: ZipFile) -> Tuple[str, str]:
    '''
    Finds and parses the VERSION file inside of an archive.

    Parameters
    ----------
    zf : ZipFile
        The zipfile object of an archive.

    Returns
    -------
    tuple of (str, str)
        The archive version and framework version of the archive.
    '''
    uuid = get_root_uuid(zf)
    version_fp = pathlib.Path(uuid) / 'VERSION'

    try:
        with zf.open(str(version_fp)) as v_fp:
            version_contents = str(v_fp.read().strip(), 'utf-8')
    except KeyError:
        raise ValueError(
            f'Malformed Archive: VERSION file for node {uuid} misplaced '
            f'or nonexistent.\nArchive {zf.filename} may be corrupt or '
            'provenance may be false.'
        )

    if not re.match(_VERSION_MATCHER, version_contents, re.MULTILINE):
        warnings.filterwarnings(
            'ignore', 'invalid escape sequence', DeprecationWarning
        )
        version_match_repr = codecs.decode(
            _VERSION_MATCHER.encode('utf-8'), 'unicode-escape'
        )
        raise ValueError(
            f'Malformed Archive: VERSION file out of spec in {zf.filename}.\n'
            f'Should match this regular expression:\n{version_match_repr}\n'
            f'Actually looks like:\n{version_contents}\n'
        )

    _, archive_version, framework_version = [
        line.strip().split()[-1] for line in
        version_contents.split(sep='\n') if line
    ]
    return (archive_version, framework_version)
