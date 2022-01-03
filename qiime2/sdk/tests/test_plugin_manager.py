# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import qiime2.plugin
import qiime2.sdk
from qiime2.plugin.plugin import SemanticTypeRecord, FormatRecord
from qiime2.sdk.plugin_manager import GetFormatFilters

from qiime2.core.testing.type import (IntSequence1, IntSequence2, IntSequence3,
                                      Mapping, FourInts, Kennel, Dog, Cat,
                                      SingleInt, C1, C2, C3, Foo, Bar, Baz,
                                      AscIntSequence, Squid, Octopus,
                                      Cuttlefish)

from qiime2.core.testing.format import (Cephalapod, IntSequenceDirectoryFormat,
                                        MappingDirectoryFormat,
                                        IntSequenceV2DirectoryFormat,
                                        IntSequenceFormatV2,
                                        IntSequenceMultiFileDirectoryFormat,
                                        FourIntsDirectoryFormat,
                                        IntSequenceFormat,
                                        RedundantSingleIntDirectoryFormat,
                                        EchoFormat,
                                        EchoDirectoryFormat,
                                        CephalapodDirectoryFormat)

from qiime2.core.testing.validator import (validator_example_null1,
                                           validate_ascending_seq,
                                           validator_example_null2)

from qiime2.core.testing.util import get_dummy_plugin


class TestPluginManager(unittest.TestCase):
    def setUp(self):
        self.plugin = get_dummy_plugin()
        # PluginManager is a singleton so there's no issue creating it again.
        self.pm = qiime2.sdk.PluginManager()

    def test_plugins(self):
        plugins = self.pm.plugins

        exp = {'dummy-plugin': self.plugin}
        self.assertEqual(plugins, exp)

    def test_validators(self):
        self.assertEqual({Kennel[Dog], Kennel[Cat], AscIntSequence, Squid,
                          Octopus, Cuttlefish},
                         set(self.pm.validators))

        self.assertEqual(
            set([r.validator for r in
                self.pm.validators[Kennel[Dog]]._validators]),
            {validator_example_null1, validator_example_null2})

        self.assertEqual(
            [r.validator for r in self.pm.validators[Kennel[Cat]]._validators],
            [validator_example_null1])

        self.assertEqual(
            [r.validator
             for r in self.pm.validators[AscIntSequence]._validators],
            [validate_ascending_seq])

    def test_type_fragments(self):
        types = self.pm.type_fragments

        exp = {
            'IntSequence1':   SemanticTypeRecord(semantic_type=IntSequence1,
                                                 plugin=self.plugin),
            'IntSequence2':   SemanticTypeRecord(semantic_type=IntSequence2,
                                                 plugin=self.plugin),
            'IntSequence3':   SemanticTypeRecord(semantic_type=IntSequence3,
                                                 plugin=self.plugin),
            'Mapping':        SemanticTypeRecord(semantic_type=Mapping,
                                                 plugin=self.plugin),
            'FourInts':       SemanticTypeRecord(semantic_type=FourInts,
                                                 plugin=self.plugin),
            'Kennel':         SemanticTypeRecord(semantic_type=Kennel,
                                                 plugin=self.plugin),
            'Dog':            SemanticTypeRecord(semantic_type=Dog,
                                                 plugin=self.plugin),
            'Cat':            SemanticTypeRecord(semantic_type=Cat,
                                                 plugin=self.plugin),
            'SingleInt':      SemanticTypeRecord(semantic_type=SingleInt,
                                                 plugin=self.plugin),
            'C1':             SemanticTypeRecord(semantic_type=C1,
                                                 plugin=self.plugin),
            'C2':             SemanticTypeRecord(semantic_type=C2,
                                                 plugin=self.plugin),
            'C3':             SemanticTypeRecord(semantic_type=C3,
                                                 plugin=self.plugin),
            'Foo':            SemanticTypeRecord(semantic_type=Foo,
                                                 plugin=self.plugin),
            'Bar':            SemanticTypeRecord(semantic_type=Bar,
                                                 plugin=self.plugin),
            'Baz':            SemanticTypeRecord(semantic_type=Baz,
                                                 plugin=self.plugin),
            'AscIntSequence': SemanticTypeRecord(semantic_type=AscIntSequence,
                                                 plugin=self.plugin),
            'Squid':          SemanticTypeRecord(semantic_type=Squid,
                                                 plugin=self.plugin),
            'Octopus':        SemanticTypeRecord(semantic_type=Octopus,
                                                 plugin=self.plugin),
            'Cuttlefish':     SemanticTypeRecord(semantic_type=Cuttlefish,
                                                 plugin=self.plugin),
        }

        self.assertEqual(types, exp)

    def test_get_semantic_types(self):
        types = self.pm.get_semantic_types()

        exp = {
            'IntSequence1': SemanticTypeRecord(semantic_type=IntSequence1,
                                               plugin=self.plugin),
            'IntSequence2': SemanticTypeRecord(semantic_type=IntSequence2,
                                               plugin=self.plugin),
            'Mapping':      SemanticTypeRecord(semantic_type=Mapping,
                                               plugin=self.plugin),
            'FourInts':     SemanticTypeRecord(semantic_type=FourInts,
                                               plugin=self.plugin),
            'Kennel[Dog]':  SemanticTypeRecord(semantic_type=Kennel[Dog],
                                               plugin=self.plugin),
            'Kennel[Cat]':  SemanticTypeRecord(semantic_type=Kennel[Cat],
                                               plugin=self.plugin),
            'SingleInt':    SemanticTypeRecord(semantic_type=SingleInt,
                                               plugin=self.plugin),
        }

        self.assertLessEqual(exp.keys(), types.keys())
        self.assertNotIn(Cat, types)
        self.assertNotIn(Dog, types)
        self.assertNotIn(Kennel, types)

    # TODO: add tests for type/directory/transformer registrations
    def test_get_formats_no_type_or_filter(self):
        exp = {
            'IntSequenceFormat':
                FormatRecord(format=IntSequenceFormat,
                             plugin=self.plugin),
            'IntSequenceDirectoryFormat':
                FormatRecord(format=IntSequenceDirectoryFormat,
                             plugin=self.plugin),
            'IntSequenceFormatV2':
                FormatRecord(format=IntSequenceFormatV2,
                             plugin=self.plugin),
            'IntSequenceV2DirectoryFormat':
                FormatRecord(format=IntSequenceV2DirectoryFormat,
                             plugin=self.plugin),
            'IntSequenceMultiFileDirectoryFormat':
                FormatRecord(format=IntSequenceMultiFileDirectoryFormat,
                             plugin=self.plugin),
            'RedundantSingleIntDirectoryFormat':
                FormatRecord(format=RedundantSingleIntDirectoryFormat,
                             plugin=self.plugin),
            'FourIntsDirectoryFormat':
                FormatRecord(format=FourIntsDirectoryFormat,
                             plugin=self.plugin),
            'EchoFormat':
                FormatRecord(format=EchoFormat,
                             plugin=self.plugin),
            'EchoDirectoryFormat':
                FormatRecord(format=EchoDirectoryFormat,
                             plugin=self.plugin),
            'MappingDirectoryFormat':
                FormatRecord(format=MappingDirectoryFormat,
                             plugin=self.plugin),
            'Cephalapod':
                FormatRecord(format=Cephalapod, plugin=self.plugin),
            'CephalapodDirectoryFormat':
                FormatRecord(format=CephalapodDirectoryFormat,
                             plugin=self.plugin),
        }

        obs = self.pm.get_formats()

        self.assertEqual(obs, exp)

    def test_get_formats_SFDF(self):
        exp = {
            'IntSequenceFormat':
                FormatRecord(format=IntSequenceFormat,
                             plugin=self.plugin),
            'IntSequenceFormatV2':
                FormatRecord(format=IntSequenceFormatV2,
                             plugin=self.plugin),
            'IntSequenceDirectoryFormat':
                FormatRecord(format=IntSequenceDirectoryFormat,
                             plugin=self.plugin),
            'IntSequenceV2DirectoryFormat':
                FormatRecord(format=IntSequenceV2DirectoryFormat,
                             plugin=self.plugin),
            'IntSequenceMultiFileDirectoryFormat':
                FormatRecord(format=IntSequenceMultiFileDirectoryFormat,
                             plugin=self.plugin)
        }

        obs = self.pm.get_formats(semantic_type='IntSequence1')

        self.assertEqual(exp, obs)

    def test_get_formats_SFDF_EXPORTABLE(self):
        exp = {
            'IntSequenceFormat':
                FormatRecord(format=IntSequenceFormat,
                             plugin=self.plugin),
            'IntSequenceFormatV2':
                FormatRecord(format=IntSequenceFormatV2,
                             plugin=self.plugin),
            'IntSequenceDirectoryFormat':
                FormatRecord(format=IntSequenceDirectoryFormat,
                             plugin=self.plugin),
            'IntSequenceV2DirectoryFormat':
                FormatRecord(format=IntSequenceV2DirectoryFormat,
                             plugin=self.plugin)
        }

        obs = self.pm.get_formats(filter=GetFormatFilters.EXPORTABLE,
                                  semantic_type=IntSequence1)

        self.assertEqual(exp, obs)

    def test_get_formats_SFDF_IMPORTABLE(self):
        exp = {
            'IntSequenceFormat':
                FormatRecord(format=IntSequenceFormat,
                             plugin=self.plugin),
            'IntSequenceDirectoryFormat':
                FormatRecord(format=IntSequenceDirectoryFormat,
                             plugin=self.plugin),
            'IntSequenceMultiFileDirectoryFormat':
                FormatRecord(format=IntSequenceMultiFileDirectoryFormat,
                             plugin=self.plugin)
        }

        obs = self.pm.get_formats(filter=GetFormatFilters.IMPORTABLE,
                                  semantic_type=IntSequence1)

        self.assertEqual(exp, obs)

    def test_get_formats_DF(self):
        exp = {
            'IntSequenceFormat':
                FormatRecord(format=IntSequenceFormat,
                             plugin=self.plugin),
            'IntSequenceFormatV2':
                FormatRecord(format=IntSequenceFormatV2,
                             plugin=self.plugin),
            'IntSequenceDirectoryFormat':
                FormatRecord(format=IntSequenceDirectoryFormat,
                             plugin=self.plugin),
            'IntSequenceV2DirectoryFormat':
                FormatRecord(format=IntSequenceV2DirectoryFormat,
                             plugin=self.plugin),
            'IntSequenceMultiFileDirectoryFormat':
                FormatRecord(format=IntSequenceMultiFileDirectoryFormat,
                             plugin=self.plugin)
        }

        obs = self.pm.get_formats(semantic_type='IntSequence3')

        self.assertEqual(exp, obs)

    def test_get_formats_DF_EXPORTABLE(self):
        exp = {
            'IntSequenceFormat':
                FormatRecord(format=IntSequenceFormat,
                             plugin=self.plugin),
            'IntSequenceDirectoryFormat':
                FormatRecord(format=IntSequenceDirectoryFormat,
                             plugin=self.plugin),
            'IntSequenceMultiFileDirectoryFormat':
                FormatRecord(format=IntSequenceMultiFileDirectoryFormat,
                             plugin=self.plugin)
        }

        obs = self.pm.get_formats(filter=GetFormatFilters.EXPORTABLE,
                                  semantic_type=IntSequence3)

        self.assertEqual(exp, obs)

    def test_get_formats_DF_IMPORTABLE(self):
        exp = {
            'IntSequenceFormatV2':
                FormatRecord(format=IntSequenceFormatV2,
                             plugin=self.plugin),
            'IntSequenceV2DirectoryFormat':
                FormatRecord(format=IntSequenceV2DirectoryFormat,
                             plugin=self.plugin),
            'IntSequenceMultiFileDirectoryFormat':
                FormatRecord(format=IntSequenceMultiFileDirectoryFormat,
                             plugin=self.plugin)
        }

        obs = self.pm.get_formats(filter=GetFormatFilters.IMPORTABLE,
                                  semantic_type=IntSequence3)

        self.assertEqual(exp, obs)

    def test_get_formats_invalid_type(self):
        with self.assertRaisesRegex(ValueError, "No formats associated"):
            self.pm.get_formats(semantic_type='Random[Frequency]')

    def test_get_formats_invalid_filter(self):
        with self.assertRaisesRegex(ValueError, "filter.*is not valid"):
            self.pm.get_formats(filter="EXPORTABLE")


if __name__ == '__main__':
    unittest.main()
