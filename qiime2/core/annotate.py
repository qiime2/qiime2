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


class Annotation():
    """General base class for all Annotation sub-classes.

    Parameters
    ----------
    name : str\n
        Name of the annotation.
        For each Result object, all Annotations must have a unique name.
        e.g. The same named Annotation can be attached to multiple Results,
        but each Result cannot contain multiple Annotations with the same name.

    Properties
    ----------
    uuid\n
        The minted uuid4 for each new Annotation that's added.
        This will be the name of each new Annotation's sub-directory within the
        `annotations` directory and is separate from any Result's uuid.

    created_at\n
        The minted date/time when an Annotation is created.
        Note that this is separate from when an Annotation is attached to a
        Results object, as this can occur at multiple times.

    Returns
    -------
    obj\n
        An instantiated Annotation of the specified sub-class.
        Note that instantiation of the Annotation base class is not supported.

    See Also
    --------
    Note

    """

    @classmethod
    def load(cls, filepath):
        """Load an Annotation.

        Parameters
        ----------
        filepath : str\n
            Path to load the Annotation from.

        Returns
        -------
        obj\n
            The instantiated Annotation sub-class.

        Raises
        ------
        ValueError\n
            If no `note.txt` file is found under the
            corresponding annotation directory.

        """
        with open(os.path.join(filepath, 'metadata.yaml'), 'r') as fh:
            meta_yaml = yaml.safe_load(fh)
            annotation_type = meta_yaml['type']

            if annotation_type == 'Note':
                annotation = Note.__new__(Note)
                # Now attach Note attrs from metadata.yaml
                annotation.name = meta_yaml['name']
                annotation.created_at = meta_yaml['created_at']
                annotation.annotation_type = meta_yaml['type']

                # Validate that `note.txt` exists
                note_fp = os.path.join(filepath, 'note.txt')
                if not os.path.exists(note_fp):
                    raise ValueError(
                        'Unable to load malformed Note with name: '
                        f'"{annotation.name}" due to missing `note.txt` file.'
                    )
                # Attach contents to Note
                else:
                    with open(note_fp, 'r') as fh:
                        annotation.contents = fh.read()

        return annotation

    # We never expect this to be hit as the base class for Annotations
    # shouldn't ever be instantiated - only the supported sub-classes.
    @property
    def annotation_type(cls):
        raise NotImplementedError

    def __init__(self, name):
        """
        Construction for an initialized Annotation.

        Attributes
        ----------
        name\n
            The user-provided name of the Annotation.

        uuid\n
            The uuid4 ID associated with the Annotation.

        created_at\n
            The datetime when the Annotation was created.

        Raises
        ------
        TypeError\n
            If the Annotation base class is instantiated.

        """
        if type(self) is Annotation:
            raise TypeError('Annotation is an abstract class'
                            ' and cannot be instantiated directly.')
        self.name = name
        self.uuid = _uuid.uuid4()
        self.created_at = datetime.now()

    def validate_name(self, name):
        """Validates that the given name is a valid Python idenitifier with the
        exception that `-` is allowed.

        Parameters
        ----------
        name : str\n
            The name to validate.

        Raises
        ------
        ValueError\n
            If the name passed in is not a valid Python identifier.
        """
        validate_name = name.replace('-', '_')
        if not validate_name.isidentifier():
            raise ValueError(f'Name "{name}" is not a valid Python identifier.'
                             ' Keys may contain `-` characters but must'
                             ' otherwise be valid Python identifiers. Python'
                             ' identifier rules may be found here'
                             ' https://www.askpython.com/python/'
                             'python-identifiers-rules-best-practices')

    def write_meta_yaml(self, annotations_dir,
                        root_result_uuid, referenced_result_uuid):
        """Write the contents of `metadata.yaml` for a given Annotation.

        Parameters
        ----------
        annotations_dir\n
            The path to the `annotations` directory within a Result object.
            Located under `provenance`.

        root_result_uuid\n
            The uuid of the Result object where an Annotation is being added.

        referenced_result_uuid\n
            The uuid of the Result object that an Annotation is referring to.
            Note that in 7.0, `root_result_uuid` and `referenced_result_uuid`
            are the same (i.e. Annotations can only refer to the same Result
            they are being attached to) but separate root and referenced uuids
            will be supported in future versions.

        Returns
        -------
        str\n
            The filepath where the Annotation's uuid-specific subdirectory
            containing the `metadata.yaml` file was written to.

        """
        # create the annotation directory if it's not already present
        # we don't care if it contains anything, just whether or not it exists
        if not os.path.exists(annotations_dir):
            os.mkdir(annotations_dir)

        # create the unique dir for a particular annotation
        annotation_uuid_dirname = \
            os.path.join(annotations_dir, str(self.uuid))
        os.mkdir(annotation_uuid_dirname)

        metadata = {
            'name': self.name,
            'id': str(self.uuid),
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
    """Note sub-class, inherits from Annotations.

    Parameters
    ----------
    text : str\n
        Inline text that will be written inside the Note's `note.txt` file.
        This parameter is optional, but either `text` OR `filepath` must be
        provided.

    filepath : str\n
        Path to a file whose contents should be written inside the Note's
        `note.txt` file.
        This parameter is optional, but either `text` OR `filepath` must be
        provided.

    Properties
    ----------
    type : Note\n
        The type of Annotation being instantiated.

    Returns
    -------
    Note : obj\n
        The instantiated Note.

    See Also
    --------
    Annotation

    """
    annotation_type = 'Note'

    # NOTE: in future versions, name will become optional & the default value
    # will be the annotation's UUID (if name isn't provided by the user)
    def __init__(self, name, *, text=None, filepath=None):
        self.validate_name(name)
        # We only want text OR filepath to be provided
        # so ensure exactly one of these gets called
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

    def write(self, annotations_dir, root_result_uuid, referenced_result_uuid):
        """Write the contents of an instantiated Note.

        Parameters
        ----------
        annotations_dir\n
            The path to the `annotations` directory within a Result object.
            Located under `provenance`.

        root_result_uuid\n
            The uuid of the Result object where an Annotation is being added.

        referenced_result_uuid\n
            The uuid of the Result object that an Annotation is referring to.
            Note that in 7.0, `root_result_uuid` and `referenced_result_uuid`
            are the same (i.e. Annotations can only refer to the same Result
            they are being attached to) but separate root and referenced uuids
            will be supported in future versions.

        See Also
        --------
        write_meta_yaml

        """
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
