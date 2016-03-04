# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import os
import uuid
import tarfile
import tempfile


# TODO use real type when it is ready
class DummyType:
    # data is simply some `bytes`

    @classmethod
    def load(cls, data_reader):
        with data_reader.get_file('data.dat', binary=True) as fh:
            return fh.read()

    @classmethod
    def save(cls, data, data_writer):
        with data_writer.create_file('data.dat', binary=True) as fh:
            fh.write(data)


# TODO check file extension on load/save and warn if not .qtf
# TODO make sure path separator code is correct for Unix and Windows systems
class Artifact:
    _readme_path = 'README.md'
    _metadata_path = 'metadata.txt'
    _data_dir = 'data'

    @classmethod
    def _get_root_dir(cls, tarfilepath):
        return os.path.splitext(os.path.basename(tarfilepath))[0]

    @classmethod
    def save(cls, data, type_, provenance, tarfilepath):
        root_dir = cls._get_root_dir(tarfilepath)

        # TODO support compression?
        with tarfile.open(tarfilepath, mode='w') as tar:
            metadata_writer = ArtifactDataWriter()
            cls._save_readme(metadata_writer)
            cls._save_metadata(type_, provenance, metadata_writer)
            metadata_writer._save_(tar, root_dir)

            data_writer = ArtifactDataWriter()
            type_.save(data, data_writer)
            data_writer._save_(tar, os.path.join(root_dir, cls._data_dir))

    @classmethod
    def _save_readme(cls, data_writer):
        with data_writer.create_file(cls._readme_path) as fh:
            fh.write(_README_TEXT)

    @classmethod
    def _save_metadata(cls, type_, provenance, data_writer):
        with data_writer.create_file(cls._metadata_path) as fh:
            fh.write(_METADATA_HEADER_TEXT)
            fh.write(cls._formatted_type(type_))
            fh.write(cls._formatted_uuid())
            fh.write(cls._formatted_provenance(provenance))

    @classmethod
    def _formatted_type(cls, type_):
        return '%s\n' % type_.__name__

    @classmethod
    def _formatted_uuid(cls):
        return '%s\n' % uuid.uuid4()

    @classmethod
    def _formatted_provenance(cls, provenance):
        return '%s\n' % provenance

    def __init__(self, tarfilepath):
        data, type_, provenance, uuid_ = self._load(tarfilepath)
        self._data = data
        self._type = type_
        self._provenance = provenance
        self._uuid = uuid_

    @property
    def data(self):
        return self._data

    @property
    def type(self):
        return self._type

    @property
    def provenance(self):
        return self._provenance

    @property
    def uuid(self):
        return self._uuid

    def _load(self, tarfilepath):
        if not tarfile.is_tarfile(tarfilepath):
            raise tarfile.ReadError(
                "%r is not a readable tar archive file" % tarfilepath)

        root_dir = self._get_root_dir(tarfilepath)
        with tarfile.open(tarfilepath, mode='r') as tar:
            metadata_reader = ArtifactDataReader(tar, root_dir)
            type_, uuid_, provenance = self._load_metadata(metadata_reader)

            # TODO lazy load data
            data_reader = ArtifactDataReader(
                tar, os.path.join(root_dir, self._data_dir))
            data = type_.load(data_reader)
            return data, type_, provenance, uuid_

    def _load_metadata(self, data_reader):
        with data_reader.get_file(self._metadata_path) as fh:
            # skip header line
            next(fh)
            type_ = self._parse_type(next(fh))
            uuid_ = self._parse_uuid(next(fh))
            provenance = self._parse_provenance(next(fh))
        return type_, uuid_, provenance

    def _parse_type(self, line):
        type_name = line.strip()

        if type_name == 'DummyType':
            return DummyType
        else:
            raise TypeError("Unrecognized artifact type %r" % type_name)

    def _parse_uuid(self, line):
        return uuid.UUID(hex=line.strip())

    def _parse_provenance(self, line):
        return line.strip()


# TODO track files, close filehandles in a special method Artifact calls
class ArtifactDataReader:
    def __init__(self, tar, data_dir):
        self._tar = tar
        self._data_dir = data_dir

    def get_paths(self):
        """Return all paths to artifact data."""
        paths = []
        for path in self._tar.getnames():
            if os.path.dirname(path) == self._data_dir:
                paths.append(os.path.basename(path))
        return paths

    def get_file(self, path, binary=False):
        try:
            filehandle = self._tar.extractfile(
                os.path.join(self._data_dir, path))
        except KeyError:
            raise FileNotFoundError("Filepath %r does not exist" % path)

        if filehandle is None:
            raise FileNotFoundError(
                "Filepath %r is not a file" % path)

        if not binary:
            filehandle = io.TextIOWrapper(filehandle)

        return filehandle


class ArtifactDataWriter:
    def __init__(self):
        # TODO pass temp dir specified by user in config
        self._tempdir = tempfile.TemporaryDirectory(
            prefix='qiime2-temp-artifact-data-')
        self._tracked_files = {}

    def create_file(self, path, binary=False):
        if self._tempdir is None:
            raise ValueError("`ArtifactDataWriter` has already been finalized")

        if path in self._tracked_files:
            raise FileExistsError("Filepath %r already exists" % path)

        if binary:
            mode = 'wb'
        else:
            mode = 'w'

        filehandle = open(os.path.join(self._tempdir.name, path), mode=mode)
        self._tracked_files[path] = filehandle
        return filehandle

    # NOTE This method is intended to be called exactly once by Artifact but
    # not plugin developers.
    def _save_(self, tar, name):
        for filehandle in self._tracked_files.values():
            filehandle.close()

        # TODO use `filter` parameter to only add files we know about
        # TODO set file metadata appropriately (e.g., owner, permissions)
        tar.add(self._tempdir.name, arcname=name)

        self._tempdir.cleanup()
        self._tempdir = None
        self._tracked_files = {}


_README_TEXT = """# QIIME 2 artifact

This tar archive file stores the data and associated metadata of a serialized
QIIME 2 artifact.

**WARNING:** This is a temporary file format used for prototyping QIIME 2. Do
not rely on this format for any reason as it will not be supported past the
prototyping stage in its current form.
"""

_METADATA_HEADER_TEXT = "# QIIME 2 artifact metadata\n"
