from enum import IntEnum
import pathlib
import warnings
from typing import Optional, Tuple
from zipfile import ZipFile

from qiime2.core.util import md5sum_directory_zip, from_checksum_format
from qiime2.core.archive.archiver import ChecksumDiff

from .util import get_root_uuid, parse_version


class ValidationCode(IntEnum):
    '''
    Codes indicating the level of validation a ProvDAG has passed.

    The code that determines which ValidationCode an archive receives is by
    necessity scattered.

    INVALID:
        One or more files are known to be missing or unparseable. Occurs
        either when checksum validation fails, or when expected files are
        absent or unparseable.
    VALIDATION_OPTOUT:
        The user opted out of checksum validation. This will be overridden by
        INVALID iff a required file is missing. In this context,
        `checksums.md5` is not required. If data files, for example, have been
        manually modified, the code will remain VALIDATION_OPTOUT, but if an
        action.yaml file is missing, INVALID will result.
    PREDATES_CHECKSUMS:
        The archive format predates the creation of checksums.md5, so full
        validation is impossible. We initially assume validity. This will be
        overridden by INVALID iff an expected file is missing or unparseable.
        If data files, for example, have been manually modified, the code will
        remain PREDATES_CHECKSUMS.
    VALID:
        The archive has passed checksum validation and is "known" to be
        valid. Md5 checksums are technically falsifiable, so this is not a
        guarantee of correctness/authenticity.
    '''
    INVALID = 0
    VALIDATION_OPTOUT = 1
    PREDATES_CHECKSUMS = 2
    VALID = 3


def validate_checksums(
    zf: ZipFile
) -> Tuple[ValidationCode, Optional[ChecksumDiff]]:
    '''
    Uses diff_checksums to validate the archive's provenance, warning the user
    if checksums.md5 is missing, or if the archive is corrupt or has been
    modified.

    Parameters
    ----------
    zf : ZipFile
        The zipfile object of the archive.

    Returns
    -------
    tuple of (ValidationCode, ChecksumDiff)
        If the checksums.md5 fle isn't present set ChecksumDiff to None
        and ValidationCode to INVALID and return.
    '''
    checksum_diff: Optional[ChecksumDiff]
    provenance_is_valid = ValidationCode.VALID

    for fp in zf.namelist():
        if 'checksums.md5' in fp:
            break
    else:
        warnings.warn(
            'The checksums.md5 file is missing from the archive. '
            'Archive may be corrupt or provenance may be false.',
            UserWarning
        )
        return ValidationCode.INVALID, None

    checksum_diff = diff_checksums(zf)
    if checksum_diff != ChecksumDiff({}, {}, {}):
        root_uuid = get_root_uuid(zf)
        warnings.warn(
            f'Checksums are invalid for Archive {root_uuid}\n'
            'Archive may be corrupt or provenance may be false.\n'
            f'Files added since archive creation: {checksum_diff.added}\n'
            'Files removed since archive creation: '
            f'{checksum_diff.removed}\n'
            'Files changed since archive creation: '
            f'{checksum_diff.changed}',
            UserWarning
        )
        provenance_is_valid = ValidationCode.INVALID

    return provenance_is_valid, checksum_diff


def diff_checksums(zf: ZipFile) -> ChecksumDiff:
    '''
    Calculates checksums for all files in an archive (except checksums.md5).
    Compares these against the checksums stored in checksums.md5, returning
    a summary ChecksumDiff.

    Parameters
    ----------
    zf : ZipFile
        The zipfile object of the archive.

    Returns
    -------
    ChecksumDiff
        A tuple of three dicts, one each for added, removed, and changed
        files. Keys are filepaths. For the added and removed dicts
        values are the checksum of the added or removed file. For the changed
        dict values are a tuple of (expected checksum, observed checksum).
    '''
    archive_version, _ = parse_version(zf)
    # TODO: don't think this is ever called
    if int(archive_version) < 5:
        return ChecksumDiff({}, {}, {})

    root_dir = pathlib.Path(get_root_uuid(zf))
    checksum_fp = str(root_dir / 'checksums.md5')
    obs = md5sum_directory_zip(zf)

    exp = {}
    for line in zf.open(checksum_fp):
        fp, checksum = from_checksum_format(str(line, 'utf-8'))
        exp[fp] = checksum

    obs_fps = set(obs)
    exp_fps = set(exp)

    added = {fp: obs[fp] for fp in obs_fps - exp_fps}
    removed = {fp: exp[fp] for fp in exp_fps - obs_fps}
    changed = {
        fp: (exp[fp], obs[fp]) for fp in exp_fps & obs_fps
        if exp[fp] != obs[fp]
    }

    return ChecksumDiff(added=added, removed=removed, changed=changed)
