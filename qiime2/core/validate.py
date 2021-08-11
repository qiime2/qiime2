# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

class ValidationChain:
    def __init__(self, concrete_type):
        self._validators = []
        self.concrete_type = concrete_type
        self._is_sorted = False

    def add_validator(self, validator):
        self._validators.append(validator)
        self._is_sorted = False

    def add_validation_object(self, *others):
        for other in others:
            self._validators += other._validators
        self._is_sorted = False


    @property
    def validators(self):
        if not self._is_sorted:
            self._validators = self._sort_validators()

        return self._validators

    def _sort_validators(self):
        """does nothing right now"""
        self._is_sorted = True

    def run_validators(self, target, validate_level: str = 'min'):
        for validator in self.validators:
            validator.validator(target, view=validator.view)
        #self._get_type_validators()
        #self._sort_validators()

        #for validator in self.sorted_validators:
        #    validator.validator(view=validator.view)
