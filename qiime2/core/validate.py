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

    def __call__(self, data, level):
        from_mt = ModelType.from_view_type(type(data))

        for record in self.validators:
            to_mt = ModelType.from_view_type(record.view)
            transformation = from_mt.make_transformation(to_mt)
            data = transformation(data)
            try:
                record.validator(data=data, level=level)
            except ValidationError:
                raise
            except Exception as e:
                raise ImplementationError("An unexpected error occured when %r"
                                          " from %r attempted to validate %r"
                                          % (record.validator.__name__,
                                             record.plugin,
                                             data)) from e

    def assert_transformation_available(self, dir_fmt):
        mt = ModelType.from_view_type(dir_fmt)

        for record in self._validators:
            mt_other = ModelType.from_view_type(record.view)
            if not mt.has_transformation(mt_other):
                raise AssertionError(
                    'Could not validate %s using %r because there was no'
                    ' transformation from %r to %r' %
                    (self.concrete_type, record.validator.__name__,
                     mt._view_name, mt_other._view_name)
                )
