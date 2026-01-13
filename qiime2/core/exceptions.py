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


class RachisWarning(UserWarning):
    """
    A custom warning that will always be displayed in the CLI, whether or not
    the --verbose flag is set.
    """
    pass
