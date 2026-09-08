# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from unittest.mock import Mock, call, patch
from typing import Union
from tempfile import TemporaryDirectory

from rachis import Artifact
from rachis import sdk
from rachis.core import util
from rachis.core.transform import (
    ModelType, NodeQueue, SearchNode, TransformType, find_transformation_path
)
from rachis.core.testing.format import (
    FirstStepFormat, SecondStepFormat, ThirdStepFormat, FourthStepFormat,
    FifthStepFormat, Cephalapod, IntSequenceFormat, IntSequenceFormatV2,
    IntSequenceDirectoryFormat, IntSequenceV2DirectoryFormat
)


class TestTransitiveTransformers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with TemporaryDirectory() as tempdir:
            cls.first_format = Artifact.import_data(
                type='FirstStep', view=tempdir
            )

        cls.int_sequence = Artifact.import_data(
            type='IntSequence1', view=[1, 2, 3]
        )

    def test_upgrade_true_only_path(self):
        """
        Path exists and each hop is `upgrade=True`.
        FirstStepFormat -> SecondStepFormat -> ThirdStepFormat
        """
        transformed = self.first_format.view(ThirdStepFormat)
        self.assertEqual(type(transformed), ThirdStepFormat)

    def test_upgrade_none_allowed_at_end_of_transformation_path(self):
        """
        The transformation from `FirstStepFormat` to `FourthStepFormat` uses
        an `upgrade=None` transformer from `ThirdStepFormat` to
        `FourthStepFormat`. This is allowed because there is only one
        `upgrade=None` hop, and it occurs at one of the ends of the path.
        """
        transformed = self.first_format.view(FourthStepFormat)
        self.assertEqual(type(transformed), FourthStepFormat)

    def test_union_transitivity(self):
        """
        Path exists between FirstStepFormat and ThirdStepFormat but not between
        FirsStepFormat and Cephalapod.
        """
        view = self.first_format.view(Union[Cephalapod, ThirdStepFormat])
        self.assertEqual(type(view), ThirdStepFormat)

    def test_int_sequence_dir_to_list(self):
        """
        Tests that unwrapping a format and then transforming does not require
        upgrading the path.
        """
        view = self.int_sequence.view(list)
        self.assertEqual(type(view), list)

    def test_lossy_transformer(self):
        """
        Tests that path with lossy steps still completes.
        First -> Second -> Third -lossy-> Fifth
        """
        view = self.first_format.view(FifthStepFormat)
        self.assertEqual(type(view), FifthStepFormat)


class TestTransitiveUpgradeSpec(unittest.TestCase):
    _implicit = object()

    def _make_path(self, *upgrades):
        '''
        Creates a dummy transformation path.

        Parameters
        ----------
        *upgrades : list[None | True | False | self._implicit]
            The upgrade types in the path. `None`, `True`, and `False`
            represent registered transformers with `upgrade=None`,
            `upgrade=True`, and `upgrade=False`, respectively. A
            `self._implicit` represents an implicit (unregistered) transformer.

        Returns
        -------
        SearchNode
            A node storing the dummy transformation path in its ancestors.
        '''
        node = SearchNode(type('Start', (), {}))

        for index, upgrade in enumerate(upgrades):
            if upgrade is self._implicit:
                record = None
                transform_type = TransformType.wrap
            else:
                record = Mock(upgrade=upgrade)
                transform_type = TransformType.registered

            node = SearchNode(
                type_=type(f'Step{index}', (), {}),
                parent=node,
                record=record,
                transform_type=transform_type,
            )

        return node

    def test_upgrade_none_only_is_valid(self):
        '''
        A path composed of a single `upgrade=None` step is valid because there
        is only one and it occurs at an end (here, both ends).
        '''
        path = self._make_path(None)
        self.assertTrue(path.validate_path())

    def test_upgrade_none_first_is_valid(self):
        '''
        A path containing a single `upgrade=None` step as the first registered
        step is valid.
        '''
        path = self._make_path(None, True)
        self.assertTrue(path.validate_path())

    def test_upgrade_none_last_is_valid(self):
        '''
        A path containging a single `upgrade=None` step as the last registered
        step is valid.
        '''
        path = self._make_path(True, None)
        self.assertTrue(path.validate_path())

    def test_implicit_steps_do_not_affect_upgrade_none_position(self):
        '''
        Implicit (unregistered) steps do not determine the "ends" of the path.
        Thus, the `upgrade=Nones` here are valid because they are the first or
        last registered transformers.
        '''
        paths = [
            self._make_path(
                self._implicit, None, True, self._implicit
            ),
            self._make_path(
                self._implicit, True, None, self._implicit
            ),
        ]
        for path in paths:
            self.assertTrue(path.validate_path())

    def test_upgrade_none_in_middle_is_invalid(self):
        '''
        A path containing an `upgrade=None` step that occurs between registered
        steps is invalid.
        '''
        path = self._make_path(True, None, True)
        self.assertFalse(path.validate_path())

    def test_multiple_upgrade_none_steps_are_invalid(self):
        '''
        A path may not contain more than one `upgrade=None` step, even if they
        occur at the ends.
        '''
        path = self._make_path(None, True, None)
        self.assertFalse(path.validate_path())

    def test_search_continues_after_invalid_target_path(self):
        '''
        Shows that an invalid target path which is explored first does not
        prevent finding a valid one later on.
        '''
        class Target:
            pass

        invalid_target = self._make_path(True, None, True)
        invalid_target.type_ = Target
        valid_target = self._make_path(True, None)
        valid_target.type_ = Target

        with patch.object(
            NodeQueue, 'pop', side_effect=[invalid_target, valid_target]
        ) as pop:
            target = find_transformation_path(object, Target)

        self.assertIs(target, valid_target)
        self.assertEqual(pop.call_count, 2)


class TestTransformationRecorder(unittest.TestCase):
    def setUp(self):
        self.pm = sdk.PluginManager()

    def _record(self, from_type, to_type):
        '''
        Mocks the recorder callable in `make_transformation`.
        '''
        recorder = Mock()
        from_mt = ModelType.from_view_type(from_type)
        to_mt = ModelType.from_view_type(to_type)
        from_mt.make_transformation(to_mt, recorder=recorder)
        return recorder

    def _call(self, record, from_type, to_type):
        '''
        Builds an expected `unittest.mock.call` object.
        '''
        input_name = util.get_view_name(from_type)
        output_name = util.get_view_name(to_type)
        return call(
            record,
            input_name=input_name,
            input_record=self.pm.views.get(input_name),
            output_name=output_name,
            output_record=self.pm.views.get(output_name),
        )

    def test_records_implicit_only_path(self):
        '''
        Asserts that a single-step implicit transformation results in one
        call to the recorder, with no `TransformerRecord`.
        '''
        recorder = self._record(
            IntSequenceFormat, IntSequenceDirectoryFormat
        )

        self.assertEqual(
            recorder.call_args_list,
            [self._call(
                None, IntSequenceFormat, IntSequenceDirectoryFormat
            )],
        )

    def test_records_identity_path(self):
        '''
        Asserts that an identity transformation results in one call to the
        recorder, with no `TransformerRecord` and equivalent start/end types.
        '''
        recorder = self._record(IntSequenceFormat, IntSequenceFormat)

        self.assertEqual(
            recorder.call_args_list,
            [self._call(None, IntSequenceFormat, IntSequenceFormat)],
        )

    def test_folds_leading_implicit_step_into_registered_step(self):
        '''
        Here we transform from `IntSequenceDirectoryFormat` to `list`. An
        implicit transformation to `IntSequenceFormat` must be made. If we
        represent the starting node in the search path as S, nodes resulting
        from an implicit transformation as I, and nodes resulting from a
        registered transformation as R, then we have:

            transformation type: S-I-R
                          index: 0-1-2

        This test shows we record the transformation as starting at S (index 0)
        and ending at R (index 2), but associating the registered transformer
        from I (index 1) to R (index 2).
        '''
        recorder = self._record(IntSequenceDirectoryFormat, list)
        record = self.pm.transformers[IntSequenceFormat][list]

        self.assertEqual(
            recorder.call_args_list,
            [self._call(record, IntSequenceDirectoryFormat, list)],
        )

    def test_folds_trailing_implicit_step_into_registered_step(self):
        '''
        See `test_folds_leading_implicit_step_into_registered_step`. Here
        we have:

            transformation type: S-R-I
                          index: 0-1-2

        This test shows that we record the transformation as starting at S
        (index 0) and ending at I (index 2), but associating the transformer
        registered from S (index 0) to R (index 1).
        '''
        recorder = self._record(list, IntSequenceDirectoryFormat)
        record = self.pm.transformers[list][IntSequenceFormat]

        self.assertEqual(
            recorder.call_args_list,
            [self._call(record, list, IntSequenceDirectoryFormat)],
        )

    def test_folds_leading_and_trailing_implicit_steps(self):
        '''
        Asserts that a path with both leading and trailing implicit steps is
        recorded as one transformation between the path endpoints, but is
        associated with the only registered transformer in the chain.

            transformation type: S-I-R-I
                          index: 0-1-2-3

        This test shows that we start at S (index 0), end at I (index 3), and
        associate the transformer from I (index 1) to R (index 2).
        '''
        recorder = self._record(
            IntSequenceDirectoryFormat, IntSequenceV2DirectoryFormat
        )
        record = self.pm.transformers[
            IntSequenceFormat][IntSequenceFormatV2]

        self.assertEqual(
            recorder.call_args_list,
            [self._call(
                record,
                IntSequenceDirectoryFormat,
                IntSequenceV2DirectoryFormat,
            )],
        )

    def test_records_each_registered_step(self):
        '''
        Asserts that when transforming from `FirstStepFormat` to
        `ThirdStepFormat` a transitive, registered hop over `SecondStepFormat`
        is made, and that the resulting two expected calls to the recorder
        appear on the mock.
        '''
        recorder = self._record(FirstStepFormat, ThirdStepFormat)

        self.assertEqual(
            recorder.call_args_list,
            [
                self._call(
                    self.pm.transformers[FirstStepFormat][SecondStepFormat],
                    FirstStepFormat,
                    SecondStepFormat,
                ),
                self._call(
                    self.pm.transformers[SecondStepFormat][ThirdStepFormat],
                    SecondStepFormat,
                    ThirdStepFormat,
                ),
            ],
        )
