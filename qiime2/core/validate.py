# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

def get_validators(concrete_type, validators):
    return validators[concrete_type]

def _sort_validators(validators):
    """does nothing right now"""
    return validators

def run_validators(concrete_type, validators):
    type_validators = get_validators(concrete_type, validators)
    sorted_validators = _sort_validators(type_validators)

    for validator in sorted_validators:
        validator.validator(view=validator.view)
