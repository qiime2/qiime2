# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2 import Metadata
from .type import Kennel, Dog, Cat
from .plugin import dummy_plugin


@dummy_plugin.register_validator(Kennel[Dog | Cat])
def test_subset_or(data: dict, validate_level):
    pass


@dummy_plugin.register_validator(Kennel[Dog])
def validator_test_null2(data: Metadata, validate_level):
    print('does know everything')
