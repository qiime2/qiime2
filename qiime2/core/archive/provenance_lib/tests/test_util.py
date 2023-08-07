import pathlib
import unittest
import zipfile

from .testing_utilities import CustomAssertions, TEST_DATA
from ..util import get_root_uuid, get_nonroot_uuid, camel_to_snake


class CamelToSnakeTests(unittest.TestCase):
    def test_camel_to_snake(self):
        some_types_n_formats = [
            'Hierarchy',  # simple
            'DistanceMatrix',  # compound
            'Bowtie2Index',  # compound w numeral
            'EMPPairedEndSequences',  # acronym
            'BIOMV210DirFmt',  # acronym w numeral
            'PCoAResults',  # weird acronym
            'PairedEndFastqManifestPhred64V2',  # compound with acronym and num
            'SampleData[Sequences]',  # bracket notation
            'FeatureData[BLAST6]',  # bracket with acronym and numeral
            'SampleData[DADA2Stats]',  # bracket complex acronym
            'FeatureData[AlignedRNASequence]',  # bracket complex acronym
            'List[FeatureTable[RelativeFrequency]]',  # made-up nested example
        ]
        exp = [
            'hierarchy',
            'distance_matrix',
            'bowtie2_index',
            'emp_paired_end_sequences',
            'biomv210_dir_fmt',
            'p_co_a_results',
            'paired_end_fastq_manifest_phred64_v2',
            'sample_data_sequences',
            'feature_data_blast6',
            'sample_data_dada2_stats',
            'feature_data_aligned_rna_sequence',
            'list_feature_table_relative_frequency',  # made-up nested example
        ]
        for og_str, exp_str in zip(some_types_n_formats, exp):
            self.assertEqual(camel_to_snake(og_str), exp_str)


class GetRootUUIDTests(unittest.TestCase):
    def test_get_root_uuid(self):
        for archive_version in TEST_DATA:
            fp = TEST_DATA[archive_version]['qzv_fp']
            exp = TEST_DATA[archive_version]['uuid']
            with zipfile.ZipFile(fp) as zf:
                self.assertEqual(exp, get_root_uuid(zf))


class GetNonRootUUIDTests(unittest.TestCase):
    def test_get_nonroot_uuid(self):
        md_example = pathlib.Path(
            'arch_root/provenance/artifacts/uuid123/metadata.yaml')
        action_example = pathlib.Path(
            'arch_root/provenance/artifacts/uuid123/action/action.yaml')
        exp = 'uuid123'

        self.assertEqual(get_nonroot_uuid(md_example), exp)
        self.assertEqual(get_nonroot_uuid(action_example), exp)


class CustomAssertionsTests(CustomAssertions):
    def test_assert_re_appears_only_once(self):
        t = ("Lick an orange. It tastes like an orange.\n"
             "The strawberries taste like strawberries!\n"
             "The snozzberries taste like snozzberries!")
        self.assertREAppearsOnlyOnce(t, 'Lick an orange')
        self.assertREAppearsOnlyOnce(t, 'tastes like')
        with self.assertRaisesRegex(AssertionError, 'Regex.*match.*orange'):
            self.assertREAppearsOnlyOnce(t, 'orange')
        with self.assertRaisesRegex(AssertionError, 'Regex.*taste like'):
            self.assertREAppearsOnlyOnce(t, 'taste like')
        with self.assertRaisesRegex(AssertionError, 'Regex.*snozzberries'):
            self.assertREAppearsOnlyOnce(t, 'snozzberries')
        with self.assertRaisesRegex(AssertionError, 'Regex.*!'):
            self.assertREAppearsOnlyOnce(t, '!')
