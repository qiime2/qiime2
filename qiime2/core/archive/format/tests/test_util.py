# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
import tempfile
import os
import zipfile

import networkx as nx

from qiime2.core.testing.type import FourInts, IntSequence1
from qiime2.core.testing.util import ArchiveTestingMixin, get_dummy_plugin
import qiime2.core.archive as archive
from qiime2.core.archive.format.util import artifact_version, parse_actions
from qiime2.sdk import Artifact


class TestArtifactVersion(unittest.TestCase, ArchiveTestingMixin):
    def setUp(self):
        prefix = "qiime2-test-temp-"
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)
        self.provenance_capture = archive.ImportProvenanceCapture()

    def test_nonexistent_archive_format(self):
        with self.assertRaisesRegex(ValueError, 'Version foo not supported'):
            with artifact_version('foo'):
                pass

    def test_write_v0_archive(self):
        fp = os.path.join(self.temp_dir.name, 'artifact_v0.qza')

        with artifact_version(0):
            artifact = Artifact._from_view(FourInts, [-1, 42, 0, 43], list,
                                           self.provenance_capture)
            artifact.save(fp)

        root_dir = str(artifact.uuid)
        # There should be no provenance
        expected = {
            'VERSION',
            'metadata.yaml',
            'data/file1.txt',
            'data/file2.txt',
            'data/nested/file3.txt',
            'data/nested/file4.txt',
        }
        self.assertArchiveMembers(fp, root_dir, expected)

        with zipfile.ZipFile(fp, mode='r') as zf:
            version = zf.read(os.path.join(root_dir, 'VERSION'))
        self.assertRegex(str(version), '^.*archive: 0.*$')

    def test_write_v1_archive(self):
        fp = os.path.join(self.temp_dir.name, 'artifact_v1.qza')

        with artifact_version(1):
            artifact = Artifact._from_view(FourInts, [-1, 42, 0, 43], list,
                                           self.provenance_capture)
            artifact.save(fp)

        root_dir = str(artifact.uuid)
        expected = {
            'VERSION',
            'metadata.yaml',
            'data/file1.txt',
            'data/file2.txt',
            'data/nested/file3.txt',
            'data/nested/file4.txt',
            'provenance/metadata.yaml',
            'provenance/VERSION',
            'provenance/action/action.yaml',
        }
        self.assertArchiveMembers(fp, root_dir, expected)

        with zipfile.ZipFile(fp, mode='r') as zf:
            version = zf.read(os.path.join(root_dir, 'VERSION'))
        self.assertRegex(str(version), '^.*archive: 1.*$')


class TestParseProvenance(unittest.TestCase, ArchiveTestingMixin):
    def setUp(self):
        prefix = "qiime2-test-temp-"
        self.plugin = get_dummy_plugin()
        self.temp_dir = tempfile.TemporaryDirectory(prefix=prefix)
        self.provenance_capture = archive.ImportProvenanceCapture()

    def test_parse_artifact_v0(self):
        fp = os.path.join(self.temp_dir.name, 'artifact_v0.qza')

        with artifact_version(0):
            artifact = Artifact._from_view(FourInts, [-1, 42, 0, 43], list,
                                           self.provenance_capture)
            artifact.save(fp)

        provenance = parse_actions(artifact)
        self.assertIsInstance(provenance, nx.DiGraph)
        # Check structure
        self.assertEqual(nx.number_of_nodes(provenance), 1)
        self.assertEqual(nx.number_of_edges(provenance), 0)
        # Check node attrs
        self.assertFalse(provenance.node[artifact.uuid])  # empty dict

    def test_parse_artifact_v1(self):
        fp = os.path.join(self.temp_dir.name, 'artifact_v1.qza')

        with artifact_version(1):
            artifact = Artifact._from_view(FourInts, [-1, 42, 0, 43], list,
                                           self.provenance_capture)
            artifact.save(fp)

        provenance = parse_actions(artifact)
        self.assertIsInstance(provenance, nx.DiGraph)
        # Check structure
        self.assertEqual(nx.number_of_nodes(provenance), 1)
        self.assertEqual(nx.number_of_edges(provenance), 0)
        # Check node attrs
        self.assertEqual(set(provenance.node[artifact.uuid]),
                         {'manifest', 'format', 'type'})

    def test_parse_two_nodes(self):
        fp = os.path.join(self.temp_dir.name, 'left.qza')
        initial = Artifact._from_view(IntSequence1, [-1, 42, 0, 43], list,
                                      self.provenance_capture)
        left, _ = self.plugin.methods['split_ints'](initial)
        left.save(fp)

        provenance = parse_actions(left)
        self.assertIsInstance(provenance, nx.DiGraph)
        # Check structure
        self.assertEqual(nx.number_of_nodes(provenance), 2)
        self.assertEqual(nx.number_of_edges(provenance), 1)
        # Check directionality
        self.assertEqual(provenance.successors(initial.uuid), [left.uuid])
        self.assertEqual(provenance.predecessors(left.uuid), [initial.uuid])
        # Check edge attrs
        self.assertEqual(set(provenance[initial.uuid][left.uuid]),
                         {'parameters', 'plugin', 'action', 'type'})

    def test_parse_three_nodes(self):
        fp = os.path.join(self.temp_dir.name, 'viz.qzv')
        initial = Artifact._from_view(IntSequence1, [-1, 42, 0, 43], list,
                                      self.provenance_capture)
        left, _ = self.plugin.methods['split_ints'](initial)
        viz, = self.plugin.visualizers['most_common_viz'](left)
        viz.save(fp)

        provenance = parse_actions(viz)
        self.assertIsInstance(provenance, nx.DiGraph)
        # Check structure
        self.assertEqual(nx.number_of_nodes(provenance), 3)
        self.assertEqual(nx.number_of_edges(provenance), 2)
        # Check directionality
        self.assertEqual(provenance.successors(initial.uuid), [left.uuid])
        self.assertEqual(provenance.successors(left.uuid), [viz.uuid])
        self.assertEqual(provenance.predecessors(viz.uuid), [left.uuid])
        self.assertEqual(provenance.predecessors(left.uuid), [initial.uuid])
        # Check edge attrs
        self.assertEqual(set(provenance[initial.uuid][left.uuid]),
                         {'parameters', 'plugin', 'action', 'type'})
        self.assertEqual(set(provenance[left.uuid][viz.uuid]),
                         {'parameters', 'plugin', 'action', 'type'})
