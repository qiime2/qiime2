# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2 import Metadata
from qiime2.plugin import ValidationError
from .type import Kennel, Dog, Cat, AscIntSequence
from .plugin import dummy_plugin


@dummy_plugin.register_validator(Kennel[Dog | Cat])
def validator_example_null1(data: dict, level):
    pass


@dummy_plugin.register_validator(Kennel[Dog])
def validator_example_null2(data: Metadata, level):
    pass


@dummy_plugin.register_validator(AscIntSequence)
def validate_ascending_seq(data: list, level):
    # landmine for testing
    if data == [2021, 8, 24]:
        raise KeyError

    prev = float('-inf')
    for number in data:
        if not number > prev:
            raise ValidationError("%s is not greater than %s" % (number, prev))
