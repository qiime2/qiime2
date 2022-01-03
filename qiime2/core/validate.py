# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.core.exceptions import ValidationError, ImplementationError
from qiime2.core.transform import ModelType
from qiime2.core.util import sorted_poset


class ValidationObject:
    r"""
    Store, sort and run all semantic validators for a for a single, complete
    semantic type(a `concrete type`).


    Attributes
    ----------
    concrete_type: SemanticType
        The semantic type for which the validators are valid for.

    """
    def __init__(self, concrete_type):
        r"""
        Create a new ValidationObject to add ValidatorRecords to.

        Parameters
        ----------
        concrete_type: semantic type
            The single, complete semantic type that the validators are to be
            associated with.

        """
        # Private Attributes
        # ------------------
        # _validators: list
        #     A list of ValidatorRecords

        # _is_sorted: Bool
        #     Tracks whether or not `_validators` has been sorted or not.

        self._validators = []
        self.concrete_type = concrete_type
        self._is_sorted = False

    def add_validator(self, validator_record):
        r"""
        Adds new validator record to plugin.

        Parameters
        ----------
        validator_record: ValidatorRecord
            ValidatorRecord is a collections.namedtuple found in
            `qiime2/plugin/plugin.py`.

        Notes
        -----
        Used by Plugin to add a `ValidatorRecord` for a new validator to a
        plugin.  Usually called through the `register_validator` decorator.

        """
        self._validators.append(validator_record)
        self._is_sorted = False

    def add_validation_object(self, *others):
        r"""
        Incorporates another validation object of the same concrete type.

        Parameters
        ----------
        *others: Any number of validation objects of the same concrete type.

        Notes
        -----
        Used to combine validation objects from different plugins. This is
        done non-heirarchically by `PluginManager` by creating a new, blank
        object for each `concrete_type` that it encounters, then adds the
        objects from each plugin.

        """
        for other in others:
            if self.concrete_type != other.concrete_type:
                raise TypeError('Unable to add ValidationObject of'
                                ' `concrete_type: %s to ValidationObject of'
                                ' `concrete_type: %s`' % (other.concrete_type,
                                                          self.concrete_type))

            self._validators += other._validators
        self._is_sorted = False

    @property
    def validators(self) -> list:
        r"""
        Public access method for the validators stored in ValidationObject.

        Returns
        -------
        list
            A sorted list of validator records.

        """
        if not self._is_sorted:
            self._sort_validators()

        return self._validators

    def _sort_validators(self):
        r"""
        Sorts validators

        Notes
        -----
        A partial order sort of the validators. The runtime for this sort is
        :math:`\theta(n^2)`. This is not a concern, as the number of
        validators present for any particular type is expected to remain
        trivially low. The validators are sorted from general to specific.

        """
        self._validators = sorted_poset(
            iterable=self._validators,
            key=lambda record: record.context,
            reverse=True)

        self._is_sorted = True

    def __call__(self, data, level):
        r"""
        Validates that provided data meets the conditions of a semantic type.

        Parameters
        ----------
        data: A view of the data to be validated.

        level: {'min', 'max'}
            specifies the level validation occurs at.

        Notes
        -----
        Use of `level` is required but the behaviour is defined in the
        individual validators.

        """
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

    def assert_transformation_available(self, data):
        r"""
        Checks that required transformations exist.

        Parameters
        ----------
        data: view
            view type of input data.

        Raises
        ------
        AssertionError
            If no transformation exists from the data view to the view
            expected by a particular validator.

        Notes
        -----
        Called by `qiime2.sdk.PluginManager._consistency_check` to ensure
        the transformers required to run the validators are defined.

        """
        mt = ModelType.from_view_type(data)

        for record in self._validators:
            mt_other = ModelType.from_view_type(record.view)
            if not mt.has_transformation(mt_other):
                raise AssertionError(
                    'Could not validate %s using %r because there was no'
                    ' transformation from %r to %r' %
                    (self.concrete_type, record.validator.__name__,
                     mt._view_name, mt_other._view_name)
                )
