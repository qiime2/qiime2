# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

class SemanticValidation:
    def __init__(self, concrete_type):
        self._validators = super.validators
        self._concrete_type = concrete_type

    def get_validators(self, concrete_type, validators):
        return validators[concrete_type]

    def _sort_validators(self, validators):
        """does nothing right now"""
        return validators

    def run_validators(self, concrete_type, validators):
        type_validators = self.get_validators(concrete_type, validators)
        sorted_validators = self._sort_validators(type_validators)

        for validator in sorted_validators:
            validator.validator(view=validator.view)
