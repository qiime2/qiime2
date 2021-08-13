# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

class ValidationObject:
    def __init__(self, concrete_type):
        self._validators = []
        self.concrete_type = concrete_type
        self._is_sorted = False

    def add_validator(self, validator_record):
        self._validators.append(validator_record)
        self._is_sorted = False

    def add_validation_object(self, *others):
        for other in others:
            self._validators += other._validators
        self._is_sorted = False


    @property
    def validators(self) -> list:
        if not self._is_sorted:
            self._validators = self._sort_validators()

        return self._validators

    def _sort_validators(self):
        """does nothing right now"""
        self._validators = self._validators
        self._is_sorted = True

    def run_validators(self, data, validate_level: str = 'min'):
        for validator in self.validators:
            validator.validator(data=data, validate_level=validate_level)

    def __call__(self, data):
        for validator in self.validators:
            validator.validator(data=data)
