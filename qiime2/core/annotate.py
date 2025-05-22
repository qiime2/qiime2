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

from collections import OrderedDict
from datetime import datetime


class Annotation():
    """General base class for all Annotation sub-classes.

    Parameters
    ----------
    name : str
        Name of the annotation.
        For each Result object, all Annotations must have a unique name.
        e.g. The same named Annotation can be attached to multiple Results,
        but each Result cannot contain multiple Annotations with the same name.

    Properties
    ----------
    uuid
        The minted uuid4 for each new Annotation that's added.
        This will be the name of each new Annotation's sub-directory within the
        `annotations` directory and is separate from any Result's uuid.

    created_at
        The minted date/time when an Annotation is created.
        Note that this is separate from when an Annotation is attached to a
        Results object, as this can occur at multiple times.

    Returns
    -------
    obj
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
        filepath : str
            Path to load the Annotation from.

        Returns
        -------
        obj
            The instantiated Annotation sub-class.

        Raises
        ------
        ValueError
            If no `note.txt` file is found under the
            corresponding annotation directory.

        """
        with open(os.path.join(filepath, 'metadata.yaml'), 'r') as fh:
            meta_yaml = yaml.safe_load(fh)
            annotation_type = meta_yaml['type']

            if annotation_type == 'Note':
                annotation = Note.__new__(Note)
                # Now attach Note attrs from metadata.yaml
                annotation.id = meta_yaml['id']
                annotation.name = meta_yaml['name']
                annotation.annotation_type = meta_yaml['type']
                annotation.created_at = meta_yaml['created_at']

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

            else:
                annotation = UnknownAnnotation.__new__(UnknownAnnotation)
                annotation.id = meta_yaml['id']
                annotation.name = meta_yaml['name']
                annotation.annotation_type = meta_yaml['type']
                annotation.created_at = meta_yaml['created_at']

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
        name
            The user-provided name of the Annotation.

        id
            The uuid4 ID associated with the Annotation.

        created_at
            The datetime when the Annotation was created.

        Raises
        ------
        TypeError
            If the Annotation base class is instantiated.

        """
        if type(self) is Annotation:
            raise TypeError('Annotation is an abstract class'
                            ' and cannot be instantiated directly.')
        self.name = name
        self.id = _uuid.uuid4()
        self.created_at = datetime.now()

    def validate_name(self, name):
        """Validates that the given name is a valid Python idenitifier with the
        exception that `-` is allowed.

        Parameters
        ----------
        name : str
            The name to validate.

        Raises
        ------
        ValueError
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

    def _write_meta_yaml(self, annotations_dir,
                         root_result_uuid, referenced_result_uuid):
        """Write the contents of `metadata.yaml` for a given Annotation.

        Parameters
        ----------
        annotations_dir
            The path to the `annotations` directory within a Result object.
            Located under `provenance`.

        root_result_uuid
            The uuid of the Result object where an Annotation is being added.

        referenced_result_uuid
            The uuid of the Result object that an Annotation is referring to.
            Note that in 7.0, `root_result_uuid` and `referenced_result_uuid`
            are the same (i.e. Annotations can only refer to the same Result
            they are being attached to) but separate root and referenced uuids
            will be supported in future versions.

        Returns
        -------
        str
            The filepath where the Annotation's uuid-specific subdirectory
            containing the `metadata.yaml` file was written to.

        """
        # create the annotation directory if it's not already present
        # we don't care if it contains anything, just whether or not it exists
        if not os.path.exists(annotations_dir):
            os.mkdir(annotations_dir)

        # create the unique dir for a particular annotation
        annotation_uuid_dirname = \
            os.path.join(annotations_dir, str(self.id))
        os.mkdir(annotation_uuid_dirname)

        metadata = OrderedDict()
        metadata['id'] = str(self.id)
        metadata['name'] = self.name
        metadata['type'] = self.annotation_type
        metadata['created_at'] = self.created_at
        metadata['root_result_uuid'] = root_result_uuid
        metadata['referenced_result_uuid'] = referenced_result_uuid

        meta_yaml = os.path.join(annotation_uuid_dirname, 'metadata.yaml')
        with open(meta_yaml, 'w') as fh:
            fh.write(yaml.dump(metadata))

        return annotation_uuid_dirname


class UnknownAnnotation(Annotation):
    """Utility sub-class that handles loading newer Annotation types on an
    older version of QIIME 2 that supports Annotations.
    """
    def __init__(*args):
        raise NotImplementedError('`UnknownAnnotation` is an abstract class'
                                  ' used for handling Annotations associated'
                                  ' with future versions of QIIME 2.'
                                  ' It should not be instantiated directly.')

    _write_meta_yaml = __init__
    _write = __init__


class Note(Annotation):
    """Note sub-class, inherits from Annotations.

    Parameters
    ----------
    text : str
        Inline text that will be written inside the Note's `note.txt` file.
        This parameter is optional, but either `text` OR `filepath` must be
        provided.

    filepath : str
        Path to a file whose contents should be written inside the Note's
        `note.txt` file.
        This parameter is optional, but either `text` OR `filepath` must be
        provided.

    Properties
    ----------
    type : Note
        The type of Annotation being instantiated.

    Returns
    -------
    Note : obj
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
        # Ensure exactly one of text or filepath is provided
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

        if text is not None:
            self.contents = text
            self._filepath = None
        else:
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

            self._filepath = str(filepath)
            self.contents = None

        # Construct Annotation class
        super().__init__(name)

    def _write(self, annotations_dir, root_result_uuid,
               referenced_result_uuid):
        """Write the contents of an instantiated Note.

        Parameters
        ----------
        annotations_dir
            The path to the `annotations` directory within a Result object.
            Located under `provenance`.

        root_result_uuid
            The uuid of the Result object where an Annotation is being added.

        referenced_result_uuid
            The uuid of the Result object that an Annotation is referring to.
            Note that in 7.0, `root_result_uuid` and `referenced_result_uuid`
            are the same (i.e. Annotations can only refer to the same Result
            they are being attached to) but separate root and referenced uuids
            will be supported in future versions.

        See Also
        --------
        _write_meta_yaml

        """
        # call _write_meta_yaml to write the stuff that's the same
        # across both input types
        annotation_uuid_dirname = \
            self._write_meta_yaml(annotations_dir,
                                  root_result_uuid,
                                  referenced_result_uuid)

        note_path = os.path.join(annotation_uuid_dirname, 'note.txt')

        if self._filepath:
            with open(self._filepath, 'rb') as fh:
                contents = fh.read()
        else:
            contents = self.contents.encode('utf-8')

        # validation for max size and parsability
        max_size = 10 * 1024 * 1024

        if len(contents) > max_size:
            raise ValueError('Note contents exceed maximum size of 10 MiB')

        try:
            contents.decode('utf-8')
        except UnicodeDecodeError:
            raise ValueError('Note contents are not valid UTF-8')

        with open(note_path, 'wb') as fh:
            fh.write(contents)
