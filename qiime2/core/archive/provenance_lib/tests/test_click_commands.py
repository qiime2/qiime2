import bibtexparser as bp
from click.testing import CliRunner
import pathlib
import tempfile
import zipfile

from ..click_commands import citations, provenance, supplement
from .test_parse import DATA_DIR, TEST_DATA
from .testing_utilities import CustomAssertions


class ReplayTests(CustomAssertions):
    def test_provenance(self):
        in_fp = TEST_DATA['5']['qzv_fp']
        in_fn = str(in_fp)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'rendered.txt'
            out_fn = str(out_fp)
            res = CliRunner().invoke(
                cli=provenance,
                args=(f"--i-in-fp {in_fn} --o-out-fp {out_fn}"))

            self.assertEqual(res.exit_code, 0)
            self.assertTrue(out_fp.is_file())

            with open(out_fn, 'r') as fp:
                rendered = fp.read()
            self.assertIn("qiime tools import", rendered)
            self.assertIn("--type 'EMPSingleEndSequences'", rendered)
            self.assertIn("--input-path <your data here>", rendered)
            self.assertIn("--output-path emp-single-end-sequences-0.qza",
                          rendered)

            self.assertREAppearsOnlyOnce(rendered, "Replay attempts.*metadata")
            self.assertRegex(rendered,
                             'The following command.*additional metadata')

            self.assertIn('qiime demux emp-single', rendered)
            self.assertIn('qiime dada2 denoise-single', rendered)
            self.assertIn('qiime phylogeny align-to-tree-mafft', rendered)
            self.assertIn(
                'recorded_metadata/diversity_core_metrics_phylogenetic_0/',
                rendered)
            self.assertIn('qiime diversity core-metrics-phylogenetic',
                          rendered)
            self.assertIn('parameter name was not found', rendered)
            self.assertIn('--?-n-jobs 1', rendered)

    def test_provenance_python(self):
        in_fp = TEST_DATA['5']['qzv_fp']
        in_fn = str(in_fp)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'rendered.txt'
            out_fn = str(out_fp)
            res = CliRunner().invoke(
                cli=provenance,
                args=(f"--i-in-fp {in_fn} --o-out-fp {out_fn} "
                      "--p-usage-driver python3"))
            self.assertEqual(res.exit_code, 0)
            self.assertTrue(out_fp.is_file())
            with open(out_fn, 'r') as fp:
                rendered = fp.read()
                self.assertIn('from qiime2 import Artifact', rendered)
                self.assertIn('import_data', rendered)
                self.assertIn('demux_actions.emp_single', rendered)
                self.assertIn('dada2_actions.denoise_single', rendered)
                self.assertIn('phylogeny_actions.align_to_tree_maf', rendered)
                self.assertIn('diversity_actions.core_metrics_phylogenetic',
                              rendered)

    def test_provenance_recurse(self):
        """
        If the directory under test is parsed recursively, two results will
        be captured from align_to_tree_mafft_fasttree instead of one.

        This is only visible in the python driver's rendering, and most users
        will never look at the underlying ProvDAG or use, so that seems like a
        reasonable way to test.
        """
        in_fp = pathlib.Path(DATA_DIR) / 'parse_dir_test'
        in_fn = str(in_fp)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'rendered.txt'
            out_fn = str(out_fp)
            res = CliRunner().invoke(
                cli=provenance,
                args=(f"--i-in-fp {in_fn} --o-out-fp {out_fn} "
                      "--p-usage-driver python3 --p-recurse"))
            self.assertEqual(res.exit_code, 0)
            self.assertTrue(out_fp.is_file())
            with open(out_fn, 'r') as fp:
                rendered = fp.read()
        self.assertRegex(
            rendered,
            '_, _, tree_0, rooted_tree_0 = phylogeny_actions.align_to_tree_m')

    def test_provenance_use_md_without_parse(self):
        in_fp = TEST_DATA['5']['qzv_fp']
        in_fn = str(in_fp)
        out_fn = 'unused_fp'
        res = CliRunner().invoke(
            cli=provenance,
            args=(f"--i-in-fp {in_fn} --o-out-fp {out_fn} "
                  "--p-no-parse-metadata --p-use-recorded-metadata"))
        self.assertEqual(res.exit_code, 1)
        self.assertIsInstance(res.exception, ValueError)
        self.assertRegex(str(res.exception), "Metadata not parsed for replay")


class ReportCitationsTests(CustomAssertions):
    def test_citations(self):
        in_fp = TEST_DATA['5']['qzv_fp']
        in_fn = str(in_fp)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'citations.bib'
            out_fn = str(out_fp)
            res = CliRunner().invoke(
                cli=citations,
                args=(f"--i-in-fp {in_fn} --o-out-fp {out_fn}"))

            self.assertEqual(res.exit_code, 0)
            self.assertTrue(out_fp.is_file())

            exp = ['action|alignment:2018.11.0|method:mafft|0',
                   'action|alignment:2018.11.0|method:mask|0',
                   'action|diversity:2018.11.0|method:beta_phylogenetic|0',
                   'action|diversity:2018.11.0|method:beta_phylogenetic|1',
                   'action|diversity:2018.11.0|method:beta_phylogenetic|2',
                   'action|diversity:2018.11.0|method:beta_phylogenetic|3',
                   'action|diversity:2018.11.0|method:beta_phylogenetic|4',
                   'action|feature-table:2018.11.0|method:rarefy|0',
                   'action|phylogeny:2018.11.0|method:fasttree|0',
                   'framework|qiime2:2018.11.0|0',
                   'plugin|dada2:2018.11.0|0',
                   'plugin|emperor:2018.11.0|0',
                   'plugin|emperor:2018.11.0|1',
                   'view|types:2018.11.0|BIOMV210DirFmt|0',
                   ]

            with open(out_fn) as bibtex_file:
                bib_database = bp.load(bibtex_file)
                self.assertEqual(len(exp), len(bib_database.entries))

            for record in set(exp):
                self.assertIn(record, bib_database.entries_dict.keys())

    def test_citations_no_deduplicate(self):
        in_fp = TEST_DATA['5']['qzv_fp']
        in_fn = str(in_fp)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'citations.bib'
            out_fn = str(out_fp)
            res = CliRunner().invoke(
                cli=citations,
                args=(f"--i-in-fp {in_fn} --o-out-fp {out_fn} "
                      "--p-no-deduplicate"))

            self.assertEqual(res.exit_code, 0)
            self.assertTrue(out_fp.is_file())

            exp = ['action|alignment:2018.11.0|method:mafft|0',
                   'action|alignment:2018.11.0|method:mask|0',
                   'action|diversity:2018.11.0|method:beta_phylogenetic|0',
                   'action|diversity:2018.11.0|method:beta_phylogenetic|1',
                   'action|diversity:2018.11.0|method:beta_phylogenetic|2',
                   'action|diversity:2018.11.0|method:beta_phylogenetic|3',
                   'action|diversity:2018.11.0|method:beta_phylogenetic|4',
                   'action|feature-table:2018.11.0|method:rarefy|0',
                   'action|phylogeny:2018.11.0|method:fasttree|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'framework|qiime2:2018.11.0|0',
                   'plugin|dada2:2018.11.0|0',
                   'plugin|dada2:2018.11.0|0',
                   'plugin|emperor:2018.11.0|0',
                   'plugin|emperor:2018.11.0|1',
                   'view|types:2018.11.0|biom.table:Table|0',
                   'view|types:2018.11.0|biom.table:Table|0',
                   'view|types:2018.11.0|BIOMV210DirFmt|0',
                   'view|types:2018.11.0|BIOMV210DirFmt|0',
                   'view|types:2018.11.0|BIOMV210DirFmt|0',
                   'view|types:2018.11.0|BIOMV210Format|0',
                   ]

            with open(out_fn) as bibtex_file:
                bib_database = bp.load(bibtex_file)
                self.assertEqual(len(exp), len(bib_database.entries))

            for record in set(exp):
                self.assertIn(record, bib_database.entries_dict.keys())


class ReproducibilitySupplementTests(CustomAssertions):
    def test_replay_supplement(self):
        in_fp = TEST_DATA['5']['qzv_fp']
        in_fn = str(in_fp)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'supplement.zip'
            out_fn = str(out_fp)
            res = CliRunner().invoke(
                cli=supplement,
                args=(f"--i-in-fp {in_fn} --o-out-fp {out_fn}"))

            self.assertEqual(res.exit_code, 0)
            self.assertTrue(out_fp.is_file())
            self.assertTrue(zipfile.is_zipfile(out_fp))

            exp = {'python3_replay.py',
                   'cli_replay.sh',
                   'citations.bib',
                   'recorded_metadata/',
                   'recorded_metadata/demux_emp_single_0/',
                   ('recorded_metadata/diversity_core_metrics_phylogenetic_0/'
                    'metadata_0.tsv'),
                   'recorded_metadata/diversity_core_metrics_phylogenetic_0/',
                   'recorded_metadata/demux_emp_single_0/barcodes_0.tsv',
                   }

            with zipfile.ZipFile(out_fp, 'r') as myzip:
                self.assertEqual(exp, set(myzip.namelist()))

    def test_replay_supplement_no_metadata_dump(self):
        """
        Confirms that metadata dumping does not occur when user opts out
        """
        in_fp = TEST_DATA['5']['qzv_fp']
        in_fn = str(in_fp)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_fp = pathlib.Path(tmpdir) / 'supplement.zip'
            out_fn = str(out_fp)
            res = CliRunner().invoke(
                cli=supplement,
                args=(f"--i-in-fp {in_fn} --o-out-fp {out_fn} "
                      "--p-no-dump-recorded-metadata"))

            self.assertEqual(res.exit_code, 0)
            self.assertTrue(out_fp.is_file())
            self.assertTrue(zipfile.is_zipfile(out_fp))

            exp = 'recorded_metadata/'

            with zipfile.ZipFile(out_fp, 'r') as myzip:
                self.assertNotIn(exp, set(myzip.namelist()))
