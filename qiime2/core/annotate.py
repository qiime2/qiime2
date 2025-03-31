# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import datetime
import os
import pathlib
import shutil
import uuid as _uuid
import yaml


class Annotation():
    @property
    def annotation_type(self):
        raise NotImplementedError
    """
    General base class for all annotation subtypes.

    Stuff that's the same here:
    - annotation uuid: uuid4 for each new annotation that's added.
        Separate from any result's uuid.

    - timestamp: current time an annotation is added to a result.

    - root_result: where is this annotation physically being added

    - referenced_result: which result (artifact/viz) this annotation refers to.
        NOTE: 7.0 only supports self-referential annotations, but this param
        is added as a note here as we would like to utilize
        referenced_result in future minor vers.

    - type
    """

    def __init__(self, root_result):
        self.uuid = _uuid.uuid4()
        self.timestamp = datetime.now()
        self.root_result = root_result

    @classmethod
    def write_meta_yaml(cls, annotations_dir, root_result_uuid):
        # annotations_dir = \
        #     (self.root_result.archive_record.root /
        #     cls.PROVENANCE_DIR /
        #     cls.ANNOTATIONS_DIR)
        if not os.path.exists(annotations_dir):
            os.mkdir(annotations_dir)

        annotation_uuid_dirname = \
            os.path.join(annotations_dir, str(_uuid.uuid4()))
        os.mkdir(annotation_uuid_dirname)

        metadata = {
            'timestamp': datetime.now(),
            'type': cls.annotation_type,
            'root_uuid': root_result_uuid
        }

        meta_yaml = os.path.join(annotation_uuid_dirname, 'metadata.yaml')
        with open(meta_yaml, 'w') as fh:
            fh.write(yaml.dump(metadata))

        return annotation_uuid_dirname


class Note(Annotation):
    annotation_type = 'Note'
    """
    Note subclass, inherits from Annotations
    - input: either inline text or a .txt file
    """
    @classmethod
    def write(cls, contents, annotations_dir, root_result_uuid):
        # make sure we're getting either a string or filepath
        if not isinstance(contents, (str, pathlib.Path)):
            raise TypeError(
                f'Unexpected input {contents}.'
                ' Accepted inputs are either a string or a filepath.'
            )
        # now call write_meta_yaml to write the stuff that's the same
        # across both input types
        annotation_uuid_dirname = \
            cls.write_meta_yaml(annotations_dir, root_result_uuid)

        # TODO: for now we're going to ignore the potential file system issues
        # associated with writing from a user provided file, but this will need
        # to be addressed prior to release via parsing/writing validation test
        # of some sort
        note_path = os.path.join(annotation_uuid_dirname, 'note.txt')

        if isinstance(contents, str):
            with open(note_path, 'w') as fh:
                fh.write(contents)
        else:
            shutil.copy(contents, note_path)


class Citation(Annotation):
    annotation_type = 'Citation'
    """
    Citation subclass, inherits from Annotations
    Should be very similar to the Citation class
    Should probably use the CitationRecord named tuple for contents
    """


class Signature(Annotation):
    annotation_type = 'Signature'
    """
    Signature subclass, inherits from Annotations
    Need to figure out what's going on in here
    """
