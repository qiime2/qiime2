# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from enum import IntEnum
import warnings
from typing import Optional, Tuple

from qiime2.core.archive.archiver import Archiver, ChecksumDiff


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
    archiver: Archiver
) -> Tuple[ValidationCode, Optional[ChecksumDiff]]:
    '''
    Uses diff_checksums to validate the archive's provenance,
    warning the user if checksums.md5/checksums.sha512 is missing,
    or if the archive is corrupt or has been modified.

    Parameters
    ----------
    result : Archiver
        The Archiver object being validated.

    Returns
    -------
    tuple of (ValidationCode, ChecksumDiff)
        If the checksums.md5/checksums.sha512 file isn't present,
        set ChecksumDiff to None and ValidationCode to INVALID and return.

    '''
    if not hasattr(archiver, 'validate_checksums'):
        return ValidationCode.PREDATES_CHECKSUMS, ChecksumDiff({}, {}, {})

    try:
        checksum_diff = archiver.validate_checksums()
    except FileNotFoundError:
        warnings.warn(
            f'The {archiver._fmt.CHECKSUM_FILE} file is missing from '
            'the archive. Archive may be corrupt or provenance may be false.',
            UserWarning
        )
        return ValidationCode.INVALID, None

    if checksum_diff != ChecksumDiff({}, {}, {}):
        warnings.warn(
            f'Checksums are invalid for Archive {archiver.uuid}\n'
            'Archive may be corrupt or provenance may be false.\n'
            f'Files added since archive creation: {checksum_diff.added}\n'
            'Files removed since archive creation: '
            f'{checksum_diff.removed}\n'
            'Files changed since archive creation: '
            f'{checksum_diff.changed}',
            UserWarning
        )
        provenance_is_valid = ValidationCode.INVALID
    else:
        provenance_is_valid = ValidationCode.VALID

    return provenance_is_valid, checksum_diff
