# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


class ValidationError(Exception):
    pass


class ImplementationError(Exception):
    pass


class QIIME2Warning(UserWarning):
    """Custom QIIME2 warning that will always be displayed"""
    pass
