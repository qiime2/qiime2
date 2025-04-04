# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import uuid as _uuid
import yaml

from datetime import datetime


# Instances of each subclass are instantiated prior to being attached
# to a Result object rather than having them be created and attached to
# the Result object directly. this allows for more flexibility
# (i.e. if the same annotation should get added to multiple Result objects)

class Annotation():
    @property
    def annotation_type(cls):
        raise NotImplementedError
    """
    General base class for all annotation subtypes.

    Stuff that's the same here:
    - annotation uuid: uuid4 for each new annotation that's added.
        Separate from any Result's uuid.

    - created_at: datetime when an annotation is created.

    - root_result: where this annotation is physically being added.

    - referenced_result: which result (artifact/viz) this annotation refers to.
        NOTE: 7.0 only supports self-referential annotations, but this param
        is added as a note here as we would like to utilize
        referenced_result in future minor vers.
    """

    @classmethod
    def load(cls, filepath):
        with open(os.path.join(filepath, 'metadata.yaml'), 'r') as fh:
            meta_yaml = yaml.safe_load(fh)
            annotation_type = meta_yaml['type']

            if annotation_type == 'Note':
                pass
                # run validation that it's an annotation and that the stuff
                # matches what we'd expect from a Note

                annotation = Note.__new__(Note)

        return annotation

    def __init__(self, name):
        self.name = name
        self.uuid = _uuid.uuid4()
        self.created_at = datetime.now()

    def write_meta_yaml(self, annotations_dir,
                        root_result_uuid, referenced_result_uuid):
        # create the annotation directory if it's not already present
        # we don't care if it exists & is empty, just whether or not it exists
        if not os.path.exists(annotations_dir):
            os.mkdir(annotations_dir)

        # create the unique dir for a particular annotation
        annotation_uuid_dirname = \
            os.path.join(annotations_dir, str(_uuid.uuid4()))
        os.mkdir(annotation_uuid_dirname)

        metadata = {
            'name': self.name,
            'created_at': self.created_at,
            'type': self.annotation_type,
            'root_result_uuid': root_result_uuid,
            'referenced_result_uuid': referenced_result_uuid
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
    def __init__(self, name, *, text=None, filepath=None):
        # We only want text OR filepath to be provided so ensure
        # exactly one of these gets called
        if text and filepath:
            raise ValueError(
                'Cannot set both `text` and `filepath` params. '
                'Please provide either inline text or a filepath only.'
            )
        if not text and not filepath:
            raise ValueError(
                'No inputs provided to either `text` or `filepath`. '
                'Please provide either inline text or a filepath.'
            )

        # Start by setting self.contents to text - this is fine if it's None
        # but will be replaced by file contents if filepath provided
        self.contents = text

        # Validate that filepath is of the correct type and file exists
        if filepath:
            if not isinstance(filepath, (str, os.PathLike)):
                raise TypeError(
                    f'Unexpected input for `filepath`: {filepath} '
                    '`filepath` should either be a '
                    '`string` or an `os.PathLike` object.'
                )
            if not os.path.exists(filepath):
                raise ValueError(
                    f'File not found from provided `filepath`: {filepath} '
                    'Double check that the provided file exists '
                    ' and is in the expected location.'
                )

            with open(filepath, 'r') as fh:
                # if we hit this branch, self.contents is now set to
                # whatever is contained in the provided file
                self.contents = fh.read()

        # Construct Annotation class
        super().__init__(name)

    # TODO: now that a Note's creation isn't immediately tied to a result,
    # need to confirm that it can be instantiated without the root_uuid
    # since this should only be attached once the Note is added
    # to a particular Result
    def write(self, annotations_dir, root_result_uuid, referenced_result_uuid):
        # call write_meta_yaml to write the stuff that's the same
        # across both input types
        annotation_uuid_dirname = \
            self.write_meta_yaml(annotations_dir,
                                 root_result_uuid,
                                 referenced_result_uuid)

        # TODO: for now we're going to ignore the potential file system issues
        # associated with writing from a user provided file, but this will need
        # to be addressed prior to release via parsing/writing validation test
        # of some sort
        note_path = os.path.join(annotation_uuid_dirname, 'note.txt')

        with open(note_path, 'w') as fh:
            fh.write(self.contents)
