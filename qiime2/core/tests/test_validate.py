# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from qiime2.core.validate import ValidationObject
from qiime2.plugin.plugin import ValidatorRecord
from qiime2.core.testing.type import IntSequence1
from qiime2.core.testing.format import IntSequenceFormat


class TestValidationObject(unittest.TestCase):

    def setUp(self):
        self.simple_int_seq = IntSequenceFormat()

        with self.simple_int_seq.open() as fh:
            fh.write('\n'.join(map(str, range(3))))
        self.simple_int_seq.validate(level='max')

    def test_initialization(self):
        validator_object = ValidationObject(IntSequence1)

        self.assertEqual(validator_object.concrete_type, IntSequence1)

    def test_add_validator(self):

        def test_validator_method(data: list, validate_level):
            pass

        test_record = ValidatorRecord(validator=test_validator_method,
                                      view=list, plugin='this_plugin',
                                      context=IntSequence1)

        validator_object = ValidationObject(IntSequence1)

        validator_object.add_validator(test_record)

        self.assertEqual(validator_object._validators,
                         [test_record])

    def test_add_validation_object(self):
        first_VO = ValidationObject(IntSequence1)
        second_VO = ValidationObject(IntSequence1)

        def first_validator(data: list, validate_level):
            pass

        def second_validator(data: list, validate_level):
            pass

        first_record = ValidatorRecord(validator=first_validator,
                                       view=list, plugin='this_plugin',
                                       context=IntSequence1)

        second_record = ValidatorRecord(validator=second_validator,
                                        view=list, plugin='this_plugin',
                                        context=IntSequence1)

        first_VO.add_validator(first_record)

        second_VO.add_validator(second_record)

        # Allows us to demonstrate add_validation_object sets _is_sorted to
        # false
        first_VO._sort_validators()

        first_VO.add_validation_object(second_VO)

        self.assertEqual(first_VO._validators, [first_record, second_record])
        self.assertFalse(first_VO._is_sorted)

    def test_public_validators_generation(self):

        validator_object = ValidationObject(IntSequence1)

        def first_validator(data: list, validate_level):
            pass

        def second_validator(data: list, validate_level):
            pass

        first_record = ValidatorRecord(validator=first_validator,
                                       view=list, plugin='this_plugin',
                                       context=IntSequence1)

        second_record = ValidatorRecord(validator=second_validator,
                                        view=list, plugin='this_plugin',
                                        context=IntSequence1)

        validator_object.add_validator(first_record)
        validator_object.add_validator(second_record)

        self.assertEqual(validator_object.validators,
                         [first_record, second_record])
        self.assertTrue(validator_object._is_sorted)

    def test_run_validators(self):

        validator_object = ValidationObject(IntSequence1)

        has_run = False

        def test_validator_method(data: list, validate_level):
            nonlocal has_run
            has_run = True
            self.assertEqual(data, [0, 1, 2])
            self.assertEqual(validate_level, 'max')

        test_record = ValidatorRecord(validator=test_validator_method,
                                      view=list, plugin='this_plugin',
                                      context=IntSequence1)

        validator_object.add_validator(test_record)

        validator_object(self.simple_int_seq, validate_level='max')

        self.assertTrue(has_run)
