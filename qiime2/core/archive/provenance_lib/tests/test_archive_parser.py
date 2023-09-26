from datetime import timedelta
import os
import networkx as nx
import pathlib
from unittest.mock import MagicMock
import pandas as pd
import tempfile
import unittest
import zipfile

import pytest

from .._checksum_validator import ChecksumDiff, ValidationCode
from .testing_utilities import (
    TestArtifacts, is_root_provnode_data, write_zip_archive
)
from ..archive_parser import (
    ProvNode, Config, _Action, _Citations, _ResultMetadata, ParserResults,
    ArchiveParser, ParserV0, ParserV1, ParserV2, ParserV3, ParserV4, ParserV5,
    ParserV6,
)
from ...provenance import MetadataInfo

from qiime2.core.testing.util import ReallyEqualMixin


class ParserVxTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_parse_root_md(self):
        for artifact in self.tas.all_artifact_versions:
            fp = artifact.filepath
            uuid = artifact.uuid
            parser = ArchiveParser.get_parser(fp)
            with zipfile.ZipFile(fp) as zf:
                root_md = parser._parse_root_md(zf, uuid)
                self.assertEqual(root_md.uuid, uuid)
                if artifact == self.tas.table_v0:
                    self.assertEqual(root_md.format, 'BIOMV210DirFmt')
                    self.assertEqual(root_md.type, 'FeatureTable[Frequency]')
                else:
                    self.assertEqual(root_md.format,
                                     'IntSequenceDirectoryFormat')
                    self.assertEqual(root_md.type, 'IntSequence1')

    def test_parse_root_md_no_md_yaml(self):
        for artifact in self.tas.all_artifact_versions:
            parser = ArchiveParser.get_parser(artifact.filepath)

            with tempfile.TemporaryDirectory() as tempdir:
                with zipfile.ZipFile(artifact.filepath) as zf:
                    zf.extractall(tempdir)

                metadata_path = os.path.join(tempdir, artifact.uuid,
                                             'metadata.yaml')
                os.remove(metadata_path)
                fn = os.path.basename(artifact.filepath)
                fp = os.path.join(tempdir, fn)
                write_zip_archive(fp, tempdir)

                with zipfile.ZipFile(fp) as zf:
                    with self.assertRaisesRegex(
                        ValueError,
                        'Malformed.*metadata'
                    ):
                        parser._parse_root_md(zf, artifact.uuid)

    @pytest.mark.filterwarnings('ignore::UserWarning')
    def test_populate_archive(self):
        for artifact in self.tas.all_artifact_versions:
            parser = ArchiveParser.get_parser(artifact.filepath)
            fp = artifact.filepath
            uuid = artifact.uuid
            version = artifact.archive_version

            if version == 0:
                with self.assertWarnsRegex(
                    UserWarning,
                    'Artifact .*prior to provenance'
                ):
                    res = parser.parse_prov(Config(), fp)

            else:
                res = parser.parse_prov(Config(), fp)
                self.assertIsInstance(res, ParserResults)
                pa_uuids = res.parsed_artifact_uuids
                self.assertIsInstance(pa_uuids, set)
                self.assertIsInstance(next(iter(pa_uuids)), str)
                self.assertIsInstance(res.prov_digraph,
                                      (type(None), nx.DiGraph))
                self.assertIsInstance(res.provenance_is_valid, ValidationCode)

                if version < 5:
                    self.assertIsInstance(res.checksum_diff, type(None))
                else:
                    self.assertIsInstance(res.checksum_diff, ChecksumDiff)

                self.assertIn(uuid, res.prov_digraph)
                self.assertIsInstance(
                    res.prov_digraph.nodes[uuid]['node_data'], ProvNode)

    def test_validate_checksums(self):
        for artifact in self.tas.all_artifact_versions:
            parser = ArchiveParser.get_parser(artifact.filepath)
            with zipfile.ZipFile(artifact.filepath) as zf:
                is_valid, diff = parser._validate_checksums(zf)
                if artifact.archive_version < 5:
                    self.assertEqual(is_valid,
                                     ValidationCode.PREDATES_CHECKSUMS)
                    self.assertEqual(diff, None)
                else:
                    self.assertEqual(is_valid, ValidationCode.VALID)
                    self.assertEqual(diff, ChecksumDiff({}, {}, {}))

    @pytest.mark.filterwarnings('ignore::UserWarning')
    def test_correct_validate_checksums_method_called(self):
        '''
        We want to confirm that parse_prov uses the local _validate_checksums
        even when it calls super().parse_prov() internally
        '''
        for artifact in self.tas.all_artifact_versions:
            parser = ArchiveParser.get_parser(artifact.filepath)
            if artifact.archive_version < 5:
                parser._validate_checksums = MagicMock(
                    # return values only here to facilitate normal execution
                    return_value=(ValidationCode.PREDATES_CHECKSUMS, None)
                )
                parser.parse_prov(Config(), artifact.filepath)
                parser._validate_checksums.assert_called_once()
            else:
                parser._validate_checksums = MagicMock(
                    return_value=(
                        ValidationCode.VALID,
                        ChecksumDiff({}, {}, {})
                    )
                )
                parser.parse_prov(Config(), artifact.filepath)
                parser._validate_checksums.assert_called_once()


class ArchiveParserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_get_parser(self):
        parsers = [
            ParserV0, ParserV1, ParserV2, ParserV3, ParserV4, ParserV5,
            ParserV6
        ]
        for artifact, parser_version in zip(
            self.tas.all_artifact_versions, parsers
        ):
            parser = ArchiveParser.get_parser(artifact.filepath)
            self.assertEqual(type(parser), parser_version)

    def test_get_parser_nonexistent_fp(self):
        fn = 'not_a_filepath.qza'
        fp = os.path.join(self.tempdir, fn)
        with self.assertRaises(FileNotFoundError):
            ArchiveParser.get_parser(fp)

    @pytest.mark.filterwarnings('ignore::UserWarning')
    def test_artifact_parser_parse_prov(self):
        with self.assertRaisesRegex(NotImplementedError, "Use a subclass"):
            ArchiveParser().parse_prov(Config(), 'doesnotmatter.txt')


class ResultMetadataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir

        cls.uuid = cls.tas.concated_ints.uuid
        md_fp = f'{cls.uuid}/provenance/metadata.yaml'
        with zipfile.ZipFile(cls.tas.concated_ints.filepath) as zf:
            cls.root_md = _ResultMetadata(zf, md_fp)

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_smoke(self):
        self.assertEqual(self.root_md.uuid, self.uuid)
        self.assertEqual(self.root_md.type, 'IntSequence1')
        self.assertEqual(self.root_md.format, 'IntSequenceDirectoryFormat')

    @pytest.mark.filterwarnings('ignore::UserWarning')
    def test_repr(self):
        exp = (f'UUID:\t\t{self.uuid}\n'
               'Type:\t\tIntSequence1\n'
               'Data Format:\tIntSequenceDirectoryFormat')
        self.assertEqual(repr(self.root_md), exp)


class ActionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir

        action_path = os.path.join(cls.tas.concated_ints_v6.uuid, 'provenance',
                                   'action', 'action.yaml')
        with zipfile.ZipFile(cls.tas.concated_ints_v6.filepath) as zf:
            cls.concat_action = _Action(zf, action_path)

        action_path = os.path.join(cls.tas.single_int.uuid, 'provenance',
                                   'action', 'action.yaml')
        with zipfile.ZipFile(cls.tas.single_int.filepath) as zf:
            cls.import_action = _Action(zf, action_path)

        action_path = os.path.join(cls.tas.pipeline_viz.uuid, 'provenance',
                                   'action', 'action.yaml')
        with zipfile.ZipFile(cls.tas.pipeline_viz.filepath) as zf:
            cls.pipeline_action = _Action(zf, action_path)

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_action_id(self):
        exp = '5035a60e-6f9a-40d4-b412-48ae52255bb5'
        self.assertEqual(self.concat_action.action_id, exp)

    def test_action_type(self):
        self.assertEqual(self.concat_action.action_type, 'method')
        self.assertEqual(self.import_action.action_type, 'import')
        self.assertEqual(self.pipeline_action.action_type, 'pipeline')

    def test_runtime(self):
        exp_t = timedelta
        exp = timedelta(microseconds=6840)
        self.assertIsInstance(self.concat_action.runtime, exp_t)
        self.assertEqual(self.concat_action.runtime, exp)

    def test_runtime_str(self):
        exp = '6840 microseconds'
        self.assertEqual(self.concat_action.runtime_str, exp)

    @pytest.mark.filterwarnings('ignore::UserWarning')
    def test_action(self):
        exp = 'concatenate_ints'
        self.assertEqual(self.concat_action.action_name, exp)

    def test_plugin(self):
        exp = 'dummy_plugin'
        self.assertEqual(self.concat_action.plugin, exp)

    '''
    Import is not handled by a plugin, so the parser provides values
    for the action_name and plugin properties not present in action.yaml
    '''
    def test_action_for_import_node(self):
        exp = 'import'
        self.assertEqual(self.import_action.action_name, exp)

    def test_plugin_for_import_node(self):
        exp = 'framework'
        self.assertEqual(self.import_action.plugin, exp)

    def test_inputs(self):
        exp = {
            'ints1': '8dea2f1a-2164-4a85-9f7d-e0641b1db22b',
            'ints2': '8dea2f1a-2164-4a85-9f7d-e0641b1db22b',
            'ints3': '7727c060-5384-445d-b007-b64b41a090ee'
        }
        self.assertEqual(self.concat_action.inputs, exp)

        exp = {}
        self.assertEqual(self.import_action.inputs, exp)

    def test_parameters(self):
        exp = {
            'int1': 7,
            'int2': 100,
        }
        self.assertEqual(self.concat_action.parameters, exp)

        exp = {}
        self.assertEqual(self.import_action.parameters, exp)

    def test_output_name(self):
        exp = 'concatenated_ints'
        self.assertEqual(self.concat_action.output_name, exp)

        exp = None
        self.assertEqual(self.import_action.output_name, exp)

    # TODO: what does format mean in this dict?
    def test_format(self):
        exp = None
        self.assertEqual(self.concat_action.format, exp)
        self.assertEqual(self.import_action.format, exp)

    def test_transformers(self):
        int_seq_dir_citation = (
            'view|dummy-plugin:0.0.0-dev|IntSequenceDirectoryFormat|0'
        )
        transformer_citation = (
            'transformer|dummy-plugin:0.0.0-dev|builtins:list'
            '->IntSequenceDirectoryFormat|0'
        )
        output_citations = [transformer_citation, int_seq_dir_citation]

        exp = {
            'inputs': {
                'ints1':
                    [{
                        'from': 'IntSequenceDirectoryFormat',
                        'to': 'builtins:list',
                        'plugin': 'dummy-plugin',
                        'citations': [int_seq_dir_citation]
                    }],
                'ints2':
                    [{
                        'from': 'IntSequenceDirectoryFormat',
                        'to': 'builtins:list',
                        'plugin': 'dummy-plugin',
                        'citations': [int_seq_dir_citation]
                    }],
                'ints3':
                    [{
                        'from': 'IntSequenceV2DirectoryFormat',
                        'to': 'builtins:list',
                        'plugin': 'dummy-plugin',
                    }],
            },
            'output': [{
                'from': 'builtins:list',
                'to': 'IntSequenceDirectoryFormat',
                'plugin': 'dummy-plugin',
                'citations': output_citations
            }]
        }
        self.assertEqual(self.concat_action.transformers, exp)

    def test_repr(self):
        exp = (
            '_Action(action_id=5035a60e-6f9a-40d4-b412-48ae52255bb5, '
            'type=method, plugin=dummy_plugin, '
            'action=concatenate_ints)'
        )
        self.assertEqual(repr(self.concat_action), exp)


class CitationsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir

        cite_strs = ['cite_none', 'cite_one', 'cite_many']
        cls.bibs = [bib+'.bib' for bib in cite_strs]
        cls.zips = [
            os.path.join(cls.tas.datadir, bib+'.zip') for bib in cite_strs
        ]

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_empty_bib(self):
        with zipfile.ZipFile(self.zips[0]) as zf:
            citations = _Citations(zf, self.bibs[0])
            self.assertEqual(len(citations.citations), 0)

    @pytest.mark.filterwarnings('ignore::UserWarning')
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
        cls.tas = TestArtifacts()
        cls.tempdir = cls.tas.tempdir

        # build root nodes for all archive format versions
        cfg = Config(parse_study_metadata=True)
        cls.nodes = {}
        for artifact in cls.tas.all_artifact_versions:
            with zipfile.ZipFile(artifact.filepath) as zf:
                all_filenames = zf.namelist()
                root_md_fnames = filter(is_root_provnode_data, all_filenames)
                root_md_fps = [pathlib.Path(fp) for fp in root_md_fnames]
                cls.nodes[str(artifact.archive_version)] = \
                    ProvNode(cfg, zf, root_md_fps)

        with zipfile.ZipFile(cls.tas.concated_ints_with_md.filepath) as zf:
            root_node_id = cls.tas.concated_ints_with_md.uuid
            all_filenames = zf.namelist()
            dag = cls.tas.concated_ints_with_md.dag
            for node in dag.nodes:
                md_path = os.path.join(
                    root_node_id, 'provenance', 'artifacts', node, 'action',
                    'metadata.tsv'
                )
                if md_path in all_filenames:
                    md_node_id = node
                else:
                    non_md_node_id = node

            # build a nonroot node without study metadata
            node_fps = [
                pathlib.Path(fp) for fp in all_filenames if
                non_md_node_id in fp and
                ('metadata.yaml' in fp or 'action.yaml' in fp
                 or 'VERSION' in fp)
            ]
            cls.nonroot_non_md_node = ProvNode(cfg, zf, node_fps)

            # build a nonroot node with study metadata
            all_filenames = zf.namelist()
            node_fps = [
                pathlib.Path(fp) for fp in all_filenames if
                md_node_id in fp and
                ('metadata.yaml' in fp or 'action.yaml' in fp
                 or 'VERSION' in fp)
            ]
            cls.nonroot_md_node = ProvNode(cfg, zf, node_fps)

            # build a root node and parse study metadata files
            root_md_fnames = filter(is_root_provnode_data, all_filenames)
            root_md_fps = [pathlib.Path(fp) for fp in root_md_fnames]
            cfg = Config(parse_study_metadata=True)
            cls.root_node_parse_md = ProvNode(cfg, zf, root_md_fps)

            # build a root node and don't parse study metadata files
            cfg = Config(parse_study_metadata=False)
            cls.root_node_dont_parse_md = ProvNode(cfg, zf, root_md_fps)

        # build a node with a collection as input
        with zipfile.ZipFile(cls.tas.int_from_collection.filepath) as zf:
            all_filenames = zf.namelist()
            root_md_fnames = filter(is_root_provnode_data, all_filenames)
            root_md_fps = [pathlib.Path(fp) for fp in root_md_fnames]
            cls.input_collection_node = ProvNode(cfg, zf, root_md_fps)

        # build a node with an optional input that defaults to None
        with zipfile.ZipFile(cls.tas.int_seq_optional_input.filepath) as zf:
            all_filenames = zf.namelist()
            root_md_fnames = filter(is_root_provnode_data, all_filenames)
            root_md_fps = [pathlib.Path(fp) for fp in root_md_fnames]
            cls.optional_input_node = ProvNode(cfg, zf, root_md_fps)

    @classmethod
    def tearDownClass(cls):
        cls.tas.free()

    def test_smoke(self):
        self.assertTrue(True)
        for node_vzn in self.nodes:
            self.assertIsInstance(self.nodes[node_vzn], ProvNode)

    def test_node_properties(self):
        # hardcoded from test data
        framework_versions = {
            '0': '2.0.5',
            '1': '2017.2.0',
            '2': '2017.9.0',
            '3': '2018.2.0',
            '4': '2018.6.0',
            '5': '2018.11.0',
            '6': '2023.5.1',
        }

        for node, archive_version in zip(
            self.nodes, [str(i) for i in range(7)]
        ):
            if archive_version == '0':
                self.assertEqual(self.nodes[node].format, 'BIOMV210DirFmt')
                self.assertEqual(self.nodes[node].type,
                                 'FeatureTable[Frequency]')
                self.assertEqual(self.nodes[node].has_provenance, False)
            else:
                self.assertEqual(self.nodes[node].format,
                                 'IntSequenceDirectoryFormat')
                self.assertEqual(self.nodes[node].type, 'IntSequence1')
                if archive_version == '1':
                    self.assertEqual(self.nodes[node].has_provenance, False)
                else:
                    self.assertEqual(self.nodes[node].has_provenance, True)

            self.assertEqual(self.nodes[node].archive_version, archive_version)
            self.assertEqual(
                self.nodes[node].framework_version,
                framework_versions[archive_version]
            )

    def test_self_eq(self):
        self.assertReallyEqual(self.nodes['5'], self.nodes['5'])

    @pytest.mark.filterwarnings('ignore::UserWarning')
    def test_eq(self):
        # Mock has no matching UUID
        mock_node = MagicMock()
        self.assertNotEqual(self.nodes['5'], mock_node)

        # Mock has bad UUID
        mock_node._uuid = 'gerbil'
        self.assertReallyNotEqual(self.nodes['5'], mock_node)

        # Matching UUIDs insufficient if classes differ
        mock_node._uuid = self.tas.concated_ints_v5.uuid
        self.assertReallyNotEqual(self.nodes['5'], mock_node)
        mock_node.__class__ = ProvNode
        self.assertReallyEqual(self.nodes['5'], mock_node)

    def test_is_hashable(self):
        exp_hash = hash(self.tas.concated_ints_v5.uuid)
        self.assertReallyEqual(hash(self.nodes['5']), exp_hash)

    def test_str(self):
        for node_vzn, artifact in zip(
            self.nodes, self.tas.all_artifact_versions
        ):
            uuid = artifact.uuid
            self.assertRegex(str(self.nodes[node_vzn]),
                             f'(?s)UUID:\t\t{uuid}.*Type.*Data Format')

    def test_repr(self):
        for node_vzn, artifact in zip(
            self.nodes, self.tas.all_artifact_versions
        ):
            uuid = artifact.uuid
            self.assertRegex(repr(self.nodes[node_vzn]),
                             f'(?s)UUID:\t\t{uuid}.*Type.*Data Format')

    def test_get_metadata_from_action(self):
        find_md = self.root_node_parse_md._get_metadata_from_Action

        # create dummy hash '0', not relevant here
        md = MetadataInfo(
            ['d5b4cf78-f5e2-44e0-aa24-66b02564e9f1'], 'metadata.tsv', '0'
        )
        action_details = {
            'parameters': [{'metadata': md}]
        }
        all_md, artifacts_as_md = find_md(action_details)
        all_exp = {'metadata': 'metadata.tsv'}
        a_as_md_exp = [{
            'artifact_passed_as_metadata':
            'd5b4cf78-f5e2-44e0-aa24-66b02564e9f1'
        }]
        self.assertEqual(all_md, all_exp)
        self.assertEqual(artifacts_as_md, a_as_md_exp)

    def test_get_metadata_from_action_with_no_params(self):
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
        self.assertEqual(type(self.nonroot_md_node.metadata), dict)
        self.assertIn('metadata', self.nonroot_md_node.metadata)
        self.assertEqual(type(self.nonroot_md_node.metadata['metadata']),
                         pd.DataFrame)

    def test_metadata_not_available_in_property_w_opt_out(self):
        self.assertEqual(self.root_node_dont_parse_md.metadata, None)

    def test_metadata_is_correct(self):
        self.assertIn('metadata', self.nonroot_md_node.metadata)

        md_data = {
            'id': ['#q2:types', '0'],
            'a': ['categorical', '42'],
        }
        md_exp = pd.DataFrame(md_data, columns=md_data.keys())
        pd.testing.assert_frame_equal(
            md_exp, self.nonroot_md_node.metadata['metadata']
        )

    def test_has_no_provenance_so_no_metadata(self):
        self.assertEqual(self.nodes['0'].has_provenance, False)
        self.assertEqual(self.nodes['0'].metadata, None)

    def test_node_has_provenance_but_no_metadata(self):
        self.assertEqual(self.nonroot_non_md_node.has_provenance, True)
        self.assertEqual(self.nonroot_non_md_node.metadata, {})

    def test_parse_metadata_for_nonroot_node(self):
        self.assertEqual(self.nonroot_md_node.has_provenance, True)
        self.assertIn('metadata', self.nonroot_md_node.metadata)

    def test_parents(self):
        actual_parent_names = []
        for parent in self.nodes['5']._parents:
            actual_parent_names += parent.keys()

        self.assertIn('ints1', actual_parent_names)
        self.assertIn('ints2', actual_parent_names)
        self.assertIn('ints3', actual_parent_names)

        self.assertEqual(len(self.nodes['5']._parents), 3)
        self.assertEqual(len(actual_parent_names), 3)

    def test_parents_no_prov(self):
        no_prov_node = self.nodes['0']
        self.assertFalse(no_prov_node.has_provenance)
        self.assertEqual(no_prov_node._parents, None)

    def test_parents_with_artifact_passed_as_md(self):
        actual_parent_names = []
        for parent in self.nonroot_md_node._parents:
            actual_parent_names += parent.keys()

        self.assertIn('ints', actual_parent_names)
        self.assertIn('artifact_passed_as_metadata', actual_parent_names)

        self.assertEqual(len(self.nonroot_md_node._parents), 2)
        self.assertEqual(len(actual_parent_names), 2)

    def test_parents_for_import_node(self):
        uuid = self.tas.single_int.uuid
        with zipfile.ZipFile(self.tas.single_int.filepath) as zf:
            required_fps = ('VERSION', 'metadata.yaml', 'action.yaml')
            import_node_fps = [
                pathlib.Path(fp) for fp in zf.namelist()
                if uuid in fp
                and any(map(lambda x: x in fp, required_fps))
            ]
            import_node = ProvNode(Config(), zf, import_node_fps)

        self.assertEqual(import_node._parents, [])

    def test_parents_collection_of_inputs(self):
        parents = self.input_collection_node._parents
        self.assertIn('int1', parents[0].keys())
        self.assertIn('int2', parents[1].keys())
        self.assertEqual(len(parents), 2)

    def test_parents_optional_input(self):
        # NOTE: The None-type input is not captured
        parents = self.optional_input_node._parents
        self.assertIn('ints', parents[0].keys())
        self.assertEqual(len(parents), 1)
