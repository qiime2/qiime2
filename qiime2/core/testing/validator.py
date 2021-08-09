# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from .type import Kennel, Dog, Cat
from .plugin import dummy_plugin


@dummy_plugin.register_validator(Kennel[Dog])
def validator_test_null(view: str):
    pass


import qiime2.core.transform as transform
@dummy_plugin.register_validator(Kennel[Dog | Cat])
def test_subset_or(view: pd.DataFrame):
    pass


@dummy_plugin.register_validator(Kennel[Dog])
def validator_test_null2(view: pd.Series):
    print('does know everything')
