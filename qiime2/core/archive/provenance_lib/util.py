import codecs
import contextlib
import os
import pathlib
import re
import warnings
import zipfile

from typing import Optional, Tuple

import qiime2
from qiime2.core.archive import Archiver


# NOTE: implemented as method on _Archive class
def get_root_uuid(zf: zipfile.ZipFile) -> str:
    """
    Returns the root UUID for a QIIME 2 Archive.

    There's no particular reason we use the first filename here.  All QIIME 2
    Artifacts store their contents in a directory named with the Artifact's
    UUID, so we can get the UUID of an artifact by taking the first part of the
    filepath of any file in the zip archive.
    """
    return pathlib.Path(zf.namelist()[0]).parts[0]


def get_nonroot_uuid(fp: pathlib.Path) -> str:
    """
    For non-root provenance files, get the Result's uuid from the path
    (avoiding the root Result's UUID which is in all paths)
    """
    if fp.name == 'action.yaml':
        uuid = fp.parts[-3]
    else:
        uuid = fp.parts[-2]
    return uuid


def camel_to_snake(name: str) -> str:
    """
    There are more comprehensive and faster ways of doing this (incl compiling)
    but it handles acronyms in semantic types nicely
    e.g. EMPSingleEndSequences -> emp_single_end_sequences
    c/o https://stackoverflow.com/a/1176023/9872253
    """
    # this will frequently be called on QIIME type expressions, so drop [ and ]
    name = re.sub(r'[\[\]]', '', name)
    # camel to snake
    name = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', name).lower()


_VERSION_MATCHER = (
    r'QIIME 2\n'
    r'archive: [0-9]{1,2}$\n'
    r'framework: '
    r'(?:20[0-9]{2}|2)\.(?:[1-9][0-2]?|0)\.[0-9](?:\.dev[0-9]?)?'
    r'(?:\+[.\w]+)?\Z'
)


def parse_version(zf: zipfile.ZipFile,
                  fp: Optional[pathlib.Path] = None) -> Tuple[str, str]:
    """Parse a VERSION file - by default uses the VERSION at archive root"""
    root_uuid = get_root_uuid(zf)
    if fp is not None:
        version_fp = fp
        node_uuid = get_nonroot_uuid(fp)
    else:
        # All files in zf start with root uuid, so we'll grab it from the first
        version_fp = pathlib.Path(root_uuid) / 'VERSION'
        node_uuid = root_uuid

    try:
        with zf.open(str(version_fp)) as v_fp:
            version_contents = str(v_fp.read().strip(), 'utf-8')
    except KeyError:

        raise ValueError(
            f"Malformed Archive: VERSION file for node {node_uuid} misplaced "
            f"or nonexistent\nArchive {zf.filename} may be corrupt or "
            "provenance may be false.")

    if not re.match(_VERSION_MATCHER, version_contents, re.MULTILINE):
        warnings.filterwarnings('ignore', 'invalid escape sequence',
                                DeprecationWarning)
        _vrsn_mtch_repr = codecs.decode(_VERSION_MATCHER.encode('utf-8'),
                                        'unicode-escape')
        raise ValueError(
            f"Malformed Archive: VERSION file out of spec in {zf.filename}\n"
            f"\nShould match this RE:\n{_vrsn_mtch_repr}\n\n"
            f"Actually looks like:\n{version_contents}\n")

    _, archive_version, frmwk_vrsn = [
        line.strip().split()[-1] for line in
        version_contents.split(sep='\n') if line]
    return (archive_version, frmwk_vrsn)


@contextlib.contextmanager
def monkeypatch_archive_version(patch_version):
    try:
        og_version = Archiver.CURRENT_FORMAT_VERSION
        Archiver.CURRENT_FORMAT_VERSION = patch_version
        yield
    finally:
        Archiver.CURRENT_FORMAT_VERSION = og_version


@contextlib.contextmanager
def monkeypatch_framework_version(patch_version):
    try:
        og_version = qiime2.__version__
        qiime2.__version__ = patch_version
        yield
    finally:
        qiime2.__version__ = og_version


def write_zip_archive(zfp, unzipped_dir):
    with zipfile.ZipFile(zfp, 'w') as zf:
        for root, dirs, files in os.walk(unzipped_dir):
            for file in files:
                path = os.path.join(root, file)
                archive_name = os.path.relpath(path, start=unzipped_dir)
                zf.write(path, arcname=archive_name)
