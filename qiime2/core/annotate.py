# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pathlib
import shutil
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

    def __init__(self, root_result):
        self.uuid = _uuid.uuid4()
        self.created_at = datetime.now()
        self.root_result = root_result
        # this is set as the root_result for now, but in the future
        # we will allow for the root and referenced results to be separate
        self.referenced_result = root_result

    def write_meta_yaml(self, name, annotations_dir, root_uuid):
        # create the annotation directory if it's not already present
        # we don't care if it exists & is empty, just whether or not it exists
        if not os.path.exists(annotations_dir):
            os.mkdir(annotations_dir)

        # create the unique dir for a particular annotation
        annotation_uuid_dirname = \
            os.path.join(annotations_dir, str(_uuid.uuid4()))
        os.mkdir(annotation_uuid_dirname)

        metadata = {
            'name': name,
            'created_at': self.created_at,
            'type': self.annotation_type,
            'root_uuid': root_uuid,
            'referenced_uuid': root_uuid
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
    # TODO: now that a Note's creation isn't immediately tied to a result,
    # need to confirm that it can be instantiated without the root_uuid
    # since this should only be attached once the Note is added
    # to a particular Result
    def write(self, name, contents, annotations_dir, root_uuid):
        # make sure we're getting either a string or filepath
        if not isinstance(contents, (str, pathlib.Path)):
            raise TypeError(
                f'Unexpected input {contents}.'
                ' Accepted inputs are either a string or a filepath.'
            )
        # now call write_meta_yaml to write the stuff that's the same
        # across both input types
        annotation_uuid_dirname = \
            self.write_meta_yaml(annotations_dir, root_uuid)

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


class Signature(Annotation):
    annotation_type = 'Signature'
    """
    Signature subclass, inherits from Annotations
    """

# TODO: depending on time, may bump this to the next minor version update
# class Citation(Annotation):
#     annotation_type = 'Citation'
#     """
#     Citation subclass, inherits from Annotations
#     Should be very similar to the Citation class
#     Should probably use the CitationRecord named tuple for contents
#     """

# TODO: depending on time, may bump this to the next minor version update
# class License(Annotation):
#     annotation_type = 'License'
#     """
#     License subclass, inherits from Annotations
#     Similar in concept to Citations
#     Will use standard Licensing formatting for validation
#     """
