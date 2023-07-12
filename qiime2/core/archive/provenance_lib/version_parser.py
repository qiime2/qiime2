import codecs
import pathlib
import re
import warnings
import zipfile
from typing import Optional, Tuple

from .util import get_nonroot_uuid, get_root_uuid

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
