from datetime import timedelta
import os
import networkx as nx
import pathlib
from unittest.mock import MagicMock
import pandas as pd
import unittest
import warnings
import zipfile

from .. import _checksum_validator
from .testing_utilities import is_root_provnode_data
from .test_parse import TEST_DATA, DATA_DIR
from ..util import UUID
from .._archive_parser import (
    ProvNode, Config, _Action, _Citations, _ResultMetadata, ParserResults,
    ArchiveParser,
)

from ...provenance import MetadataInfo
from .testing_utilities import ReallyEqualMixin


class ParserVxTests(unittest.TestCase):
    def test_parse_root_md(self):
        for archive_version in TEST_DATA:
            fp = TEST_DATA[archive_version]['qzv_fp']
            root_uuid = TEST_DATA[archive_version]['uuid']
            parser = TEST_DATA[archive_version]['parser']()
            with zipfile.ZipFile(fp) as zf:
                root_md = parser._parse_root_md(zf, root_uuid)
                self.assertEqual(root_md.uuid, root_uuid)
                self.assertEqual(root_md.type,  'Visualization')
                self.assertEqual(root_md.format, None)

    def test_parse_root_md_no_md_yaml(self):
        qzv_no_root_md = os.path.join(DATA_DIR, 'no_root_md_yaml.qzv')
        for archive_version in TEST_DATA:
            root_uuid = TEST_DATA[archive_version]['uuid']
            parser = TEST_DATA[archive_version]['parser']()
            with zipfile.ZipFile(qzv_no_root_md) as zf:
                with self.assertRaisesRegex(ValueError, 'Malformed.*metadata'):
                    parser._parse_root_md(zf, root_uuid)

    def test_populate_archive(self):
        for archive_version in TEST_DATA:
            qzv_fp = TEST_DATA[archive_version]['qzv_fp']
            root_uuid = TEST_DATA[archive_version]['uuid']
            parser = TEST_DATA[archive_version]['parser']()
            if archive_version == '0':
                with self.assertWarnsRegex(
                    UserWarning,
                        'Artifact 0b8b47.*prior to provenance'):
                    res = parser.parse_prov(Config(), qzv_fp)
            else:
                res = parser.parse_prov(Config(), qzv_fp)

                self.assertIsInstance(res, ParserResults)
                pa_uuids = res.parsed_artifact_uuids
                self.assertIsInstance(pa_uuids, set)
                self.assertIsInstance(next(iter(pa_uuids)), UUID)
                self.assertIsInstance(res.prov_digraph,
                                      (type(None), nx.DiGraph))
                self.assertIsInstance(res.provenance_is_valid,
                                      _checksum_validator.ValidationCode)
                exp_diff_type = (type(None) if int(archive_version[0]) < 5
                                 else _checksum_validator.ChecksumDiff)
                self.assertIsInstance(res.checksum_diff, exp_diff_type)

                # Does this archive have the expected number of Results?
                self.assertEqual(len(res.prov_digraph),
                                 TEST_DATA[archive_version]['n_res'])
                # Is the root UUID a key a node in the DiGraph?
                self.assertIn(root_uuid, res.prov_digraph)
                # Is contents keyed on uuids, containing ProvNodes?
                self.assertIsInstance(
                    res.prov_digraph.nodes[root_uuid]['node_data'],
                    ProvNode)

    def test_validate_checksums(self):
        for archive_version in TEST_DATA:
            with zipfile.ZipFile(TEST_DATA[archive_version]['qzv_fp']) as zf:
                parser = TEST_DATA[archive_version]['parser']()
                is_valid, diff = parser._validate_checksums(zf)
                self.assertEqual(is_valid,
                                 TEST_DATA[archive_version]['prov_is_valid'])
                self.assertEqual(diff,
                                 TEST_DATA[archive_version]['checksum'])

    def test_correct_validate_checksums_method_called(self):
        # We want to confirm that parse_prov uses the local _validate_checksums
        # even when it calls super().parse_prov() internally
        for archive_version in TEST_DATA:
            parser = TEST_DATA[archive_version]['parser']()
            parser._validate_checksums = MagicMock(
                # return values only here to facilitate normal execution
                return_value=(TEST_DATA[archive_version]['prov_is_valid'],
                              TEST_DATA[archive_version]['checksum']))
            qzv_fp = TEST_DATA[archive_version]['qzv_fp']
            if archive_version == '0':
                # supress warning from parsing provenance for a v0 ProvDAG
                uuid = TEST_DATA['0']['uuid']
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        'ignore',  f'Art.*{uuid}.*prior')
                    parser.parse_prov(Config(), qzv_fp)
                    parser._validate_checksums.assert_called_once()
            else:
                parser.parse_prov(Config(), qzv_fp)
                parser._validate_checksums.assert_called_once()


class ArchiveParserTests(unittest.TestCase):
    def test_get_parser(self):
        for version in TEST_DATA:
            fp = os.path.join(DATA_DIR, TEST_DATA[version]['qzv_fp'])
            parser = ArchiveParser.get_parser(fp)
            self.assertIsInstance(parser, TEST_DATA[version]['parser'])

    def test_get_parser_nonexistent_fp(self):
        fn = 'not_a_filepath.qza'
        fp = os.path.join(DATA_DIR, fn)
        with self.assertRaisesRegex(
                FileNotFoundError, f'ArchiveParser.*{fp}'):
            ArchiveParser.get_parser(fp)

    def test_get_parser_insufficient_permissions(self):
        fn = 'not_a_zip.txt'
        fp = os.path.join(DATA_DIR, fn)
        # Hack this, because it's impossible to commit a file with 0 read perms
        os.chmod(fp, 0o000)
        with self.assertRaisesRegex(
                PermissionError, f"ArchiveParser.*denied.*{fn}"):
            ArchiveParser.get_parser(fp)
        os.chmod(fp, 0o644)

    def test_get_parser_not_a_zip_archive(self):
        fn = 'not_a_zip.txt'
        fp = os.path.join(DATA_DIR, fn)
        with self.assertRaisesRegex(
                zipfile.BadZipFile, "ArchiveParser.*File is not a zip file"):
            ArchiveParser.get_parser(fp)

    def test_artifact_parser_parse_prov(self):
        with self.assertRaisesRegex(NotImplementedError, "Use a subclass"):
            ArchiveParser().parse_prov(Config(), 'doesnotmatter.txt')


class ResultMetadataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        v5_qzv = TEST_DATA['5']['qzv_fp']
        cls.v5_uuid = TEST_DATA['5']['uuid']
        md_fp = f'{cls.v5_uuid}/provenance/metadata.yaml'
        with zipfile.ZipFile(str(v5_qzv)) as zf:
            cls.v5_root_md = _ResultMetadata(zf, md_fp)

    def test_smoke(self):
        self.assertEqual(self.v5_root_md.uuid, self.v5_uuid)
        self.assertEqual(self.v5_root_md.type, 'Visualization')
        self.assertEqual(self.v5_root_md.format, None)

    def test_repr(self):
        exp = (f'UUID:\t\t{self.v5_uuid}\n'
               'Type:\t\tVisualization\n'
               'Data Format:\tNone')
        self.assertEqual(repr(self.v5_root_md), exp)


class ActionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root_action_fp = os.path.join(DATA_DIR,
                                      'action_emperor_root_node_v5.zip')
        import_action_fp = os.path.join(DATA_DIR, 'action_import_v5.zip')
        with zipfile.ZipFile(root_action_fp) as zf:
            cls.act = _Action(zf, 'action.yaml')

        with zipfile.ZipFile(import_action_fp) as zf:
            cls.imp_act = _Action(zf, 'action.yaml')

    def test_action_id(self):
        exp = '5bc4b090-abbc-46b0-a219-346c8026f7d7'
        self.assertEqual(self.act.action_id, exp)

    def test_action_type(self):
        exp = 'pipeline'
        self.assertEqual(self.act.action_type, exp)

    def test_runtime(self):
        exp_t = timedelta
        exp = timedelta(seconds=2, microseconds=17110)
        self.assertIsInstance(self.act.runtime, exp_t)
        self.assertEqual(self.act.runtime, exp)

    def test_runtime_str(self):
        exp = '2 seconds, and 17110 microseconds'
        self.assertEqual(self.act.runtime_str, exp)

    def test_action(self):
        exp = 'core_metrics_phylogenetic'
        self.assertEqual(self.act.action_name, exp)

    def test_plugin(self):
        exp = 'diversity'
        self.assertEqual(self.act.plugin, exp)

    def test_inputs(self):
        action_exp = {'table': '706b6bce-8f19-4ae9-b8f5-21b14a814a1b',
                      'phylogeny': 'ad7e5b50-065c-4fdd-8d9b-991e92caad22'}
        import_exp = {}
        self.assertEqual(self.act.inputs, action_exp)
        self.assertEqual(self.imp_act.inputs, import_exp)

    def test_parameters(self):
        action_exp = {'sampling_depth': 1000,
                      'metadata': MetadataInfo(input_artifact_uuids=[],
                                               relative_fp='metadata.tsv'),
                      'n_jobs_or_threads': 'auto'}
        import_exp = {}
        self.assertEqual(self.act.parameters, action_exp)
        self.assertEqual(self.imp_act.parameters, import_exp)

    def test_output_name(self):
        action_exp = 'unweighted_unifrac_emperor'
        import_exp = None
        self.assertEqual(self.act.output_name, action_exp)
        self.assertEqual(self.imp_act.output_name, import_exp)

    def test_format(self):
        action_exp = None
        import_exp = 'EMPSingleEndDirFmt'
        self.assertEqual(self.act.format, action_exp)
        self.assertEqual(self.imp_act.format, import_exp)

    def test_transformers(self):
        action_exp = None
        import_exp = {'output': [{'from': 'EMPSingleEndDirFmt',
                                  'to': 'EMPSingleEndDirFmt'}]}
        self.assertEqual(self.act.transformers, action_exp)
        self.assertEqual(self.imp_act.transformers, import_exp)

    def test_repr(self):
        exp = ('_Action(action_id=5bc4b090-abbc-46b0-a219-346c8026f7d7, '
               'type=pipeline, plugin=diversity, '
               'action=core_metrics_phylogenetic)')
        self.assertEqual(repr(self.act), exp)

    # NOTE: Import is not handled by a plugin, so the parser provides values
    # for the action_name and plugin properties not present in action.yaml
    def test_action_for_import_node(self):
        exp = 'import'
        self.assertEqual(self.imp_act.action_name, exp)

    def test_plugin_for_import_node(self):
        exp = 'framework'
        self.assertEqual(self.imp_act.plugin, exp)


class CitationsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cite_strs = ['cite_none', 'cite_one', 'cite_many']
        cls.bibs = [bib+'.bib' for bib in cite_strs]
        cls.zips = [os.path.join(DATA_DIR, bib+'.zip') for bib in cite_strs]

    def test_empty_bib(self):
        with zipfile.ZipFile(self.zips[0]) as zf:
            citations = _Citations(zf, self.bibs[0])
            # Is the _citations dict empty?
            self.assertFalse(len(citations.citations))

    def test_citation(self):
        with zipfile.ZipFile(self.zips[1]) as zf:
            exp = 'framework'
            citations = _Citations(zf, self.bibs[1])
            for key in citations.citations:
                self.assertRegex(key, exp)

    def test_many_citations(self):
        exp = ['2020.6.0.dev0', 'unweighted_unifrac.+0',
               'unweighted_unifrac.+1', 'unweighted_unifrac.+2',
               'unweighted_unifrac.+3', 'unweighted_unifrac.+4',
               'BIOMV210DirFmt', 'BIOMV210Format']
        with zipfile.ZipFile(self.zips[2]) as zf:
            citations = _Citations(zf, self.bibs[2])
            for i, key in enumerate(citations.citations):
                self.assertRegex(key, exp[i])

    def test_repr(self):
        exp = ("Citations(['framework|qiime2:2020.6.0.dev0|0'])")
        with zipfile.ZipFile(self.zips[1]) as zf:
            citations = _Citations(zf, self.bibs[1])
            self.assertEqual(repr(citations), exp)


class ProvNodeTests(unittest.TestCase, ReallyEqualMixin):
    @classmethod
    def setUpClass(cls):
        cfg = Config(parse_study_metadata=True)
        # Build root nodes for all archive format versions
        cls.nodes = dict()
        for k in list(TEST_DATA):
            with zipfile.ZipFile(TEST_DATA[k]['qzv_fp']) as zf:
                all_filenames = zf.namelist()
                root_md_fnames = filter(is_root_provnode_data, all_filenames)
                root_md_fps = [pathlib.Path(fp) for fp in root_md_fnames]
                cls.nodes[k] = ProvNode(cfg, zf, root_md_fps)

        # Build a minimal node in which Artifacts are passed as metadata
        # NOTE: This file breaks some assumptions about qzas for simplicity.
        # e.g. it doesn't contain provenance data for its results passed as md
        # and may need to be replaced if these nodes adopt new behavior.
        filename = pathlib.Path('minimal_v4_artifact_as_md.zip')
        # This manufactured v4 test archive "happens" to have the same root
        # UUID as the standard v5 archive
        pfx = pathlib.Path(TEST_DATA['5']['uuid'])
        artifact_as_md_fp = os.path.join(DATA_DIR, filename)
        with zipfile.ZipFile(artifact_as_md_fp) as zf:
            cls.art_as_md_node = ProvNode(
                cfg,
                zf,
                [pfx / 'VERSION',
                 pfx / 'metadata.yaml',
                 pfx / 'provenance/action/action.yaml']
                )

        # Build a nonroot node without study metadata
        with zipfile.ZipFile(TEST_DATA['5']['qzv_fp']) as zf:
            node_id = '3b7d36ff-37ab-4ac2-958b-6a547d442bcf'
            all_filenames = zf.namelist()
            node_fps = [
                pathlib.Path(fp) for fp in all_filenames if
                node_id in fp and
                ('metadata.yaml' in fp or 'action.yaml' in fp
                 or 'VERSION' in fp
                 )]
            cls.nonroot_non_md_node = ProvNode(cfg, zf, node_fps)

            # Build a nonroot node with study metadata
            node_id = '0af08fa8-48b7-4c6a-83c6-e0f766156343'
            all_filenames = zf.namelist()
            node_fps = [
                pathlib.Path(fp) for fp in all_filenames if
                node_id in fp and
                ('metadata.yaml' in fp or 'action.yaml' in fp
                 or 'VERSION' in fp
                 )]
            cls.nonroot_md_node = ProvNode(cfg, zf, node_fps)

            # Build a root node and don't parse study metadata files
            node_id = TEST_DATA['5']['uuid']
            root_md_fnames = filter(is_root_provnode_data, zf.namelist())
            root_md_fps = [pathlib.Path(fp) for fp in root_md_fnames]
            cfg = Config(parse_study_metadata=False)
            cls.dont_parse_md_files_node = ProvNode(cfg, zf, root_md_fps)

        # build a node with a collection as input
        with zipfile.ZipFile(os.path.join(DATA_DIR, 'merged_tbls.qza')) as zf:
            all_filenames = zf.namelist()
            root_md_fnames = filter(is_root_provnode_data, all_filenames)
            root_md_fps = [pathlib.Path(fp) for fp in root_md_fnames]
            cls.input_collection_node = ProvNode(cfg, zf, root_md_fps)

        # build a node with an optional input that defaults to None
        with zipfile.ZipFile(
                os.path.join(DATA_DIR, 'optional_input_none.qzv')) as zf:
            all_filenames = zf.namelist()
            root_md_fnames = filter(is_root_provnode_data, all_filenames)
            root_md_fps = [pathlib.Path(fp) for fp in root_md_fnames]
            cls.optional_input_node = ProvNode(cfg, zf, root_md_fps)

    def test_smoke(self):
        self.assertTrue(True)
        for node_vzn in self.nodes:
            self.assertIsInstance(self.nodes[node_vzn], ProvNode)

    def test_properties_with_viz(self):
        for node in self.nodes:
            self.assertEqual(self.nodes[node]._uuid, TEST_DATA[node]['uuid'])
            self.assertEqual(self.nodes[node].type, 'Visualization')
            self.assertEqual(self.nodes[node].format, None)
            self.assertEqual(self.nodes[node].archive_version,
                             TEST_DATA[node]['av'])
            self.assertEqual(self.nodes[node].framework_version,
                             TEST_DATA[node]['fwv'])
            self.assertEqual(self.nodes[node].has_provenance,
                             TEST_DATA[node]['has_prov'])

    def test_self_eq(self):
        self.assertReallyEqual(self.nodes['5'], self.nodes['5'])

    def test_eq(self):
        # Mock has no matching UUID
        mock_node = MagicMock()
        self.assertNotEqual(self.nodes['5'], mock_node)

        # Mock has bad UUID
        mock_node._uuid = 'gerbil'
        self.assertReallyNotEqual(self.nodes['5'], mock_node)

        # Matching UUIDs insufficient if classes differ
        mock_node._uuid = TEST_DATA['5']['uuid']
        self.assertReallyNotEqual(self.nodes['5'], mock_node)
        mock_node.__class__ = ProvNode
        self.assertReallyEqual(self.nodes['5'], mock_node)

    def test_is_hashable(self):
        exp_hash = hash(TEST_DATA['5']['uuid'])
        self.assertReallyEqual(hash(self.nodes['5']), exp_hash)

    def test_str(self):
        for node_vzn in self.nodes:
            uuid = TEST_DATA[node_vzn]['uuid']
            self.assertRegex(repr(self.nodes[node_vzn]),
                             f'(?s)UUID:\t\t{uuid}.*Type.*Data Format')

    def test_repr(self):
        for node_vzn in self.nodes:
            uuid = TEST_DATA[node_vzn]['uuid']
            self.assertRegex(repr(self.nodes[node_vzn]),
                             f'(?s)UUID:\t\t{uuid}.*Type.*Data Format')

    def test_archive_version(self):
        for node_vzn in self.nodes:
            self.assertEqual(self.nodes[node_vzn].archive_version,
                             TEST_DATA[node_vzn]['av'])

    def test_framework_version(self):
        for node_vzn in self.nodes:
            self.assertEqual(self.nodes[node_vzn].framework_version,
                             TEST_DATA[node_vzn]['fwv'])

    def test_get_metadata_from_action(self):
        find_md = self.nodes['5']._get_metadata_from_Action
        md1 = MetadataInfo([], 'some_metadata.tsv')
        md2 = MetadataInfo(['301b4'], 'other_metadata.tsv')
        md3 = MetadataInfo(['4154', '5555b'], 'merged_metadata.tsv')
        action_details = \
            {'parameters':
                [
                 {'some_param': 'foo'},
                 {'arbitrary_metadata_name': md1},
                 {'other_metadata': md2},
                 {'double_md': md3},
                 ]}
        all_md, artifacts_as_md = find_md(action_details)
        all_exp = {'arbitrary_metadata_name': 'some_metadata.tsv',
                   'other_metadata': 'other_metadata.tsv',
                   'double_md': 'merged_metadata.tsv',
                   }
        a_as_md_exp = [{'artifact_passed_as_metadata': '301b4'},
                       {'artifact_passed_as_metadata': '4154'},
                       {'artifact_passed_as_metadata': '5555b'},
                       ]
        self.assertEqual(all_md, all_exp)
        self.assertEqual(artifacts_as_md, a_as_md_exp)

    def test_get_metadata_from_action_with_actual_node(self):
        find_md = self.nodes['5']._get_metadata_from_Action
        all_md, artifacts_as_md = find_md(
            self.nodes['5'].action._action_details)
        exp = {'metadata': 'metadata.tsv'}
        self.assertEqual(all_md, exp)
        self.assertEqual(artifacts_as_md, [])

    def test_get_metadata_from_action_with_no_params(self):
        # Not sure which of these test data are possible in action.yaml, but
        # we'll check them both just in case
        find_md = self.nodes['5']._get_metadata_from_Action
        action_details = \
            {'parameters': []}
        all_md, artifacts_as_md = find_md(action_details)
        self.assertEqual(all_md, {})
        self.assertEqual(artifacts_as_md, [])

        action_details = {'non-parameters-key': 'here is a thing'}
        all_md, artifacts_as_md = find_md(action_details)
        self.assertEqual(all_md, {})
        self.assertEqual(artifacts_as_md, [])

    def test_metadata_available_in_property(self):
        self.assertEqual(type(self.nodes['5'].metadata), dict)
        self.assertIn('metadata', self.nodes['5'].metadata)
        self.assertEqual(type(self.nodes['5'].metadata['metadata']),
                         pd.DataFrame)

    def test_metadata_not_available_in_property_w_opt_out(self):
        self.assertEqual(self.dont_parse_md_files_node.metadata, None)

    def test_metadata_is_correct(self):
        # Were parameter names captured correctly?
        self.assertIn('sample_metadata', self.art_as_md_node.metadata)
        self.assertIn('feature_metadata', self.art_as_md_node.metadata)

        # Does sample metadata look right?
        s_m_data = {'sampleid': ['#q2:types', 's_id_123'],
                    'barcodeSequence': ['categorical', 'TACCGCTTCTTC'],
                    'isGerbil': ['categorical', 'totally']}
        s_m_exp = pd.DataFrame(s_m_data, columns=s_m_data.keys())
        pd.testing.assert_frame_equal(
            s_m_exp,
            self.art_as_md_node.metadata['sample_metadata'])

        # Does feature metadata look right?
        f_m_data = {'Feature ID': ['#q2:types', 'feature_id_123'],
                    'Taxon': ['categorical', 'd__Bacteria; p__Firmicutes'],
                    'Confidence': ['categorical', '0.9']}
        f_m_exp = pd.DataFrame(f_m_data, columns=f_m_data.keys())
        pd.testing.assert_frame_equal(
            f_m_exp,
            self.art_as_md_node.metadata['feature_metadata'])

    def test_has_no_provenance_so_no_metadata(self):
        self.assertEqual(self.nodes['0'].has_provenance, False)
        self.assertEqual(self.nodes['0'].metadata, None)

    def test_node_has_provenance_but_no_metadata(self):
        self.assertIn('3b7d36ff', self.nonroot_non_md_node._uuid)
        self.assertEqual(self.nonroot_non_md_node.has_provenance, True)
        self.assertEqual(self.nonroot_non_md_node.metadata, {})

    def test_parse_metadata_for_nonroot_node(self):
        self.assertIn('0af08fa8', self.nonroot_md_node._uuid)
        self.assertEqual(self.nonroot_md_node.has_provenance, True)
        self.assertIn('metadata', self.nonroot_md_node.metadata)

    def test_parents(self):
        exp = [{'table': '89af91c0-033d-4e30-8ac4-f29a3b407dc1'},
               {'phylogeny': 'bce3d09b-e296-4f2b-9af4-834db6412429'}]
        self.assertEqual(self.nodes['5']._parents, exp)

    def test_parents_no_prov(self):
        no_prov_node = self.nodes['0']
        self.assertFalse(no_prov_node.has_provenance)
        self.assertEqual(no_prov_node._parents, None)

    def test_parents_with_artifact_passed_as_md(self):
        exp = [{'tree': 'e710bdc5-e875-4876-b238-5451e3e8eb46'},
               {'feature_table': 'abc22fdc-e7fa-4976-a980-8f2ff8c4bb58'},
               {'pcoa': '1ed04b10-d29c-495f-996e-3d4db89434d2'},
               {'artifact_passed_as_metadata':
                '415409a4-371d-4c69-9433-e3eaba5301b4'},
               ]
        actual = self.art_as_md_node._parents
        self.assertEqual(actual, exp)

    def test_parents_for_import_node(self):
        with zipfile.ZipFile(TEST_DATA['5']['qzv_fp']) as zf:
            import_node_id = 'a35830e1-4535-47c6-aa23-be295a57ee1c'
            reqd_fps = ('VERSION', 'metadata.yaml', 'action.yaml')
            import_node_fps = [
                pathlib.Path(fp) for fp in zf.namelist()
                if import_node_id in fp
                and any(map(lambda x: x in fp, reqd_fps))
                ]
            import_node = ProvNode(Config(), zf, import_node_fps)

        self.assertEqual(import_node._parents, [])

    def test_parents_for_table_with_collection_of_inputs(self):
        exp = [{'tables_0': "84898e39-f6e0-44bb-8fa1-6df2f330af68"},
               {'tables_1': "0be6c7be-ad84-4417-9f1c-cade0a8a9b58"}]
        parents = self.input_collection_node._parents
        self.assertEqual(parents, exp)

    def test_parents_for_table_with_optional_input(self):
        # NOTE: The None-type input is not captured
        exp = [{'artifact_passed_as_metadata':
                "1cc80a0b-1415-49b1-9d58-bd9394e5f613"}]
        parents = self.optional_input_node._parents
        self.assertEqual(parents, exp)
