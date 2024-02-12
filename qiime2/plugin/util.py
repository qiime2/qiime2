# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import psutil

from qiime2.core.transform import ModelType


def transform(data, *, from_type=None, to_type):
    from_type = type(data) if from_type is None else from_type

    from_model_type = ModelType.from_view_type(from_type)
    to_model_type = ModelType.from_view_type(to_type)
    transformation = from_model_type.make_transformation(to_model_type)

    return transformation(data)


def get_available_cores(one_less: bool = False):
    '''
    Finds the number of currently available (logical) cores. Useful for plugins
    that need to convert a 0 to a concrete number of cores when 0 is not
    supported by the underlying/called software.

    Parameters
    ----------
    one_less : bool
        Whether to return one less than the total number of cores available.

    Returns
    -------
    int
        The number of cores to be requested.
    '''
    cpus = psutil.cpu_count()
    if cpus is not None and cpus > 1:
        return cpus - one_less

    return 1
