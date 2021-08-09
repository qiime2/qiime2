# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

class SemanticValidation:
    def __init__(self, validation_target, concrete_type, validators, validate_level):
        self.data = validation_target
        self.validators = validators
        self.concrete_type = concrete_type
        self.validate_level = validate_level

    def _get_type_validators(self):
        self.type_validators = self.validators[self.concrete_type]

    def _sort_validators(self):
        """does nothing right now"""
        self.sorted_validators = self.type_validators
        pass

    def run_validators(self):
        self._get_type_validators()
        self._sort_validators()

        for validator in self.sorted_validators:
            validator.validator(view=validator.view)
