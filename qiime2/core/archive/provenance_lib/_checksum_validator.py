from enum import IntEnum
import pathlib
import warnings
import zipfile
from typing import Optional, Tuple

from qiime2.core.util import md5sum_directory_zip, from_checksum_format
from qiime2.core.archive.archiver import ChecksumDiff

from .util import get_root_uuid, parse_version


class ValidationCode(IntEnum):
    """
    Codes indicating the level of validation a ProvDAG has passed.

    The code that determines what ValidationCode an archive receives is by
    necessity scattered. Though not ideal, this is probably the best
    "central" location to keep information on when these codes will occur.

    INVALID: one or more files are known to be missing or unparseable. Occurs
        either when checksum validation fails, or when expected files are
        absent or unparseable.
    VALIDATION_OPTOUT: The user opted out of checksum validation. This will be
        overridden by INVALID iff a required file is missing. In this context,
        `checksums.md5` is not required. If data files, for example, have been
        manually modified, the code will remain VALIDATION_OPTOUT, but if an
        action.yaml file is missing, INVALID will result.
    PREDATES_CHECKSUMS: The archive format predates the creation of
        checksums.md5, so full validation is impossible. We initially assume
        validity. This will be overridden by INVALID iff an expected file is
        missing or unparseable.  If data files, for example, have been manually
        modified, the code will remain PREDATES_CHECKSUMS
    VALID: The archive has passed checksum validation and is "known" to be
        valid. Md5 checksums are technically falsifiable, so this is not a
        guarantee of correctness/authenticity. It would, however, require a
        significant and unlikely effort at falsification of results to render
        this untrue.
    """
    INVALID = 0                 # Archive is known to be invalid
    VALIDATION_OPTOUT = 1       # User opted out of validation
    PREDATES_CHECKSUMS = 2      # v0-v4 cannot be validated, so assume validity
    VALID = 3                   # Archive known to be valid


def validate_checksums(zf: zipfile.ZipFile) -> Tuple[ValidationCode,
                                                     Optional[ChecksumDiff]]:
    """
    Uses diff_checksums to validate the archive's provenance, warning the user
    if checksums.md5 is missing, or if the archive is corrupt/has been modified

    Returns a (ValidationCode, ChecksumDiff) tuple. For archive formats prior
    to v5, the ChecksumDiff will be empty b/c checksums.md5 does not exist.

    The returned ChecksumDiff will be None iff checksums.md5 should be present
    (b/c v5+) but is missing.
    """
    checksum_diff: Optional[ChecksumDiff]
    provenance_is_valid = ValidationCode.VALID

    # One broad try/except here saves us more down the call stack
    try:
        checksum_diff = diff_checksums(zf)
        if checksum_diff != ChecksumDiff({}, {}, {}):
            # self._result_md may not have been parsed yet, so get uuid
            root_uuid = get_root_uuid(zf)
            warnings.warn(
                f"Checksums are invalid for Archive {root_uuid}\n"
                "Archive may be corrupt or provenance may be false"
                ".\n"
                f"Files added since archive creation: {checksum_diff.added}\n"
                "Files removed since archive creation: "
                f"{checksum_diff.removed}\n"
                "Files changed since archive creation: "
                f"{checksum_diff.changed}", UserWarning)
            provenance_is_valid = ValidationCode.INVALID
    # zipfiles KeyError if file not found. warn if checksums.md5 is missing
    # and return ChecksumDiff=None
    except KeyError as err:
        warnings.warn(
            str(err).strip('"') +
            ". Archive may be corrupt or provenance may be false",
            UserWarning)
        provenance_is_valid = ValidationCode.INVALID
        checksum_diff = None

    return (provenance_is_valid, checksum_diff)


def diff_checksums(zf: zipfile.ZipFile) -> ChecksumDiff:
    """
    Calculates checksums for all files in an archive (excepting checksums.md5)
    Compares these against the checksums stored in checksums.md5, returning
    a summary ChecksumDiff

    For archive formats prior to v5, returns an empty ChecksumDiff b/c
    checksums.md5 does not exist

    Code adapted from qiime2/core/archive/archiver.py
    """
    archive_version, _ = parse_version(zf)
    if int(archive_version) < 5:
        return ChecksumDiff({}, {}, {})

    root_dir = pathlib.Path(get_root_uuid(zf))
    checksum_filename = root_dir / 'checksums.md5'
    obs = md5sum_directory_zip(zf)
    exp = dict(from_checksum_format(str(line, 'utf-8')) for line in
               zf.open(str(checksum_filename))
               )
    obs_keys = set(obs)
    exp_keys = set(exp)

    added = {x: obs[x] for x in obs_keys - exp_keys}
    removed = {x: exp[x] for x in exp_keys - obs_keys}
    changed = {x: (exp[x], obs[x]) for x in exp_keys & obs_keys
               if exp[x] != obs[x]}

    return ChecksumDiff(added=added, removed=removed, changed=changed)
