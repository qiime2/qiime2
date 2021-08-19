# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.core.exceptions import ValidationError, ImplementationError
from qiime2.core.transform import ModelType


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
        # self._validators = self._validators
        self._is_sorted = True
        return self._validators

    def __call__(self, data, validate_level):

        from_mt = ModelType.from_view_type(type(data))

        for validator in self.validators:
            to_mt = ModelType.from_view_type(validator.view)
            transformation = from_mt.make_transformation(to_mt)
            data = transformation(data)
            try:
                validator.validator(data=data, validate_level=validate_level)
            except ValidationError:
                raise
            except Exception as e:
                raise ImplementationError("An unexpected error occured when %s"
                                          " from %s attempted to validate %r"
                                          % (validator.validator,
                                             validator.plugin,
                                             data)) from e
