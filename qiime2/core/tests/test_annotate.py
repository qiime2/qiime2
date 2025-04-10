# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest

from qiime2.core.annotate import Note
from qiime2.core.testing.type import FourInts
from qiime2.sdk.result import Artifact, Result


class TestAnnotations(unittest.TestCase):
    def setUp(self):
        # Create Artifact
        test_dir = tempfile.TemporaryDirectory(prefix='qiime2-test-temp-')
        saved_artifact = Artifact.import_data(FourInts, [-1, 42, 0, 43])
        saved_artifact.save(os.path.join(test_dir.name, 'artifact.qza'))

        self.artifact = \
            Result.load(os.path.join(test_dir.name, 'artifact.qza'))

        # Create Notes
        self.note1 = Note(name='mynote', text='my special text')
        self.note2 = Note(name='mynote', text='my other special text')
        self.note3 = Note(name='mynote3', text='my extra special text')

    def test_add_annotation_roundtrip(self):
        # confirm that annotations starts as an empty list
        self.assertEqual(self.artifact._annotations, [])

        # add note1 to ints1 artifact
        self.artifact.add_annotation(self.note1)

        # check that there's exactly one annotation entry
        self.assertEqual(len(self.artifact._annotations), 1)

        for annotation in self.artifact.iter_annotations():
            self.assertEqual(annotation.name, 'mynote')
            self.assertEqual(annotation.annotation_type, 'Note')
            self.assertEqual(annotation.contents, 'my special text')

            # TODO: pull out annotation dir name for systematic usage
            # add test to assert contents of annotation dir

# class specs ##
# instantiating a base class annotation (E)
# incorrectly instantiating note (E)
# correctly instantiating two notes for use in downstream testing

# endpoints ##
# add_annotation
# - attempt to add another annotation with non-unique name (E)

# get_annotation
# - attempt to get an annotation from a result wo any annotations (E)
# - get annotation from above result w/added annotation
# - attempt to get an annotation w/name not found (E)

# iter_annotations
# - get annotations from result with multiple notes & check
# expected details match per annotation
# - attempt to call on result wo any annotations (E)

# remove_annotation
# - attempt to remove annotation from result wo any annotations (E)
# - attempt to remove annotation wo found name (E)
# - remove annotation from above result
