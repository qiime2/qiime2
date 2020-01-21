# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.core.transform import ModelType


def transform(data, *, from_type=None, to_type):
    from_type = type(data) if from_type is None else from_type

    from_model_type = ModelType.from_view_type(from_type)
    to_model_type = ModelType.from_view_type(to_type)
    transformation = from_model_type.make_transformation(to_model_type)

    return transformation(data)
