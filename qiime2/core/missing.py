# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import numpy as np

from .enan import make_nan_with_payload as _make_nan_with_payload
from .enan import get_payload_from_nan as _get_payload_from_nan


def _encode_terms(namespace, enum):
    namespace = _NAMESPACE_LOOKUP.index(namespace)

    def encode(x):
        try:
            code = enum.index(x)
        except ValueError:
            return x
        return _make_nan_with_payload(code, namespace=namespace)

    return encode


def _insdc_missing(series):
    return series.apply(
        _encode_terms('INSDC:missing', (
            'not applicable', 'missing', 'not collected', 'not provided',
            'restricted access'))
    )


def _q2_omitted(series):
    return series


def _q2_error(series):
    if series.isna().any():
        raise ValueError("NOT ALLOWED")
    return series


BUILTIN_MISSING = {
    'INSDC:missing': _insdc_missing,
    'q2:omitted': _q2_omitted,
    'q2:error': _q2_error
}


# list index reflects the term used, the "q2:" enums don't apply here, since
# they aren't actually encoded in the NaNs
_NAMESPACE_LOOKUP = ['INSDC:missing']
DEFAULT_MISSING = 'q2:omitted'


def series_encode_missing(series: pd.Series, enumeration: str) -> pd.Series:
    if type(enumeration) is not str:
        TypeError("Wrong type for `enumeration`, expected string")
    try:
        encoder = BUILTIN_MISSING[enumeration]
    except KeyError:
        raise ValueError("Unknown enumeration: %r, (available: %r)"
                         % (enumeration, list(BUILTIN_MISSING.keys())))

    new = encoder(series)
    if series.dtype == object and new.isna().all():
        # return to categorical of all missing values
        return new.astype(object)
    return new


def series_extract_missing(series: pd.Series) -> pd.Series:
    def _decode(x):
        if np.issubdtype(type(x), np.floating) and np.isnan(x):
            code, namespace = _get_payload_from_nan(x)
            if namespace is None:
                return x
            elif namespace == 255:
                raise ValueError("Custom enumerations are not yet supported")
            else:
                try:
                    _, enum = BUILTIN_MISSING[_NAMESPACE_LOOKUP[namespace]]
                except (IndexError, KeyError):
                    return x

            try:
                return enum[code]
            except IndexError:
                return x

        return x

    missing = series[series.isna()]
    return missing.apply(_decode)
