import unittest
from unittest.mock import Mock, call
from typing import Union
from tempfile import TemporaryDirectory

from rachis import Artifact
from rachis import sdk
from rachis.core import util
from rachis.core.transform import ModelType
from rachis.core.testing.format import (
    FirstStepFormat, SecondStepFormat, ThirdStepFormat, FourthStepFormat,
    FifthStepFormat, Cephalapod, IntSequenceFormat, IntSequenceFormatV2,
    IntSequenceDirectoryFormat, IntSequenceV2DirectoryFormat)


class TestTransitiveTransfomrers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with TemporaryDirectory() as tempdir:
            cls.first_format = Artifact.import_data(
                type='FirstStep', view=tempdir
            )

        cls.int_sequence = Artifact.import_data(
            type='IntSequence1', view=[1, 2, 3]
        )

    def test_first_to_third(self):
        """
        Path exists and is upgraded.
        FirstStepFormat -> SecondStepFormat -> ThirdStepFormat
        """
        view = self.first_format.view(ThirdStepFormat)
        self.assertEqual(type(view), ThirdStepFormat)

    def test_first_to_fourth_fails(self):
        """
        Path exists but is not upgraded.
        FirstStepFormat -> SecondStepFormat -> ThirdStepFormat -None-> Fourth
        """
        with self.assertRaisesRegex(Exception, 'No transformation from'):
            self.first_format.view(FourthStepFormat)

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
