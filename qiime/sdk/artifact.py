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
import yaml
import importlib
import time
import datetime

from qiime.sdk.type import Type


class Artifact:
    _readme_path = 'README.md'
    _metadata_path = 'metadata.yaml'
    _data_dir = 'data'

    @classmethod
    def _get_root_dir(cls, tarfilepath):
        return os.path.splitext(os.path.basename(tarfilepath))[0]

    @classmethod
    def save(cls, data, type_, provenance, tarfilepath):
        root_dir = cls._get_root_dir(tarfilepath)

        with tarfile.open(tarfilepath, mode='w') as tar:
            metadata_writer = ArtifactDataWriter()
            cls._save_readme(metadata_writer)
            cls._save_metadata(type_, provenance, metadata_writer)
            metadata_writer._save_(tar, root_dir)

            data_writer = ArtifactDataWriter()
            type_().save(data, data_writer)
            data_writer._save_(tar, os.path.join(root_dir, cls._data_dir))

    @classmethod
    def _save_readme(cls, data_writer):
        with data_writer.create_file(cls._readme_path) as fh:
            fh.write(_README_TEXT)

    @classmethod
    def _save_metadata(cls, type_, provenance, data_writer):
        if not type_().is_concrete():
            raise TypeError("%r is not a concrete type. Only concrete types "
                            "can be saved." % type_)
        type_exp = repr(type_)
        # TODO collapse imports with common prefix
        imports = [':'.join([path, name]) for name, path
                   in type_().get_imports()]

        with data_writer.create_file(cls._metadata_path) as fh:
            yaml.safe_dump({
                'UUID': cls._formatted_uuid(),
                'provenance': cls._formatted_provenance(provenance),
                'imports': imports,
                'type': type_exp
            }, stream=fh)

    @classmethod
    def _formatted_uuid(cls):
        return str(uuid.uuid4())

    @classmethod
    def _formatted_provenance(cls, provenance):
        if provenance is None:
            return ('No provenance available because this artifact was '
                    'generated independently of QIIME.')
        else:
            return {'workflow-UUID': str(provenance.workflow_uuid),
                    'artifact-uuids': provenance.artifact_uuids,
                    # TODO: need to store parameter types (or is it enough
                    # that the workflow template reference could get us that?)
                    'parameters': provenance.parameters,
                    'runtime': cls._formatted_runtime(provenance),
                    'workflow-template-reference': provenance.workflow_template_reference}

    @classmethod
    def _formatted_runtime(cls, provenance):
        stop_time = time.time()
        run_time = stop_time - provenance.start_time
        datetime_start = datetime.datetime.fromtimestamp(provenance.start_time)
        datetime_stop = datetime.datetime.fromtimestamp(stop_time)
        return {'run-time-seconds': run_time,
                'execution-start-time': datetime_start.isoformat(),
                'execution-stop-time': datetime_stop.isoformat()}

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

            data_reader = ArtifactDataReader(
                tar, os.path.join(root_dir, self._data_dir))
            data = type_().load(data_reader)
            return data, type_, provenance, uuid_

    def _load_metadata(self, data_reader):
        with data_reader.get_file(self._metadata_path) as fh:
            metadata = yaml.safe_load(fh)
            uuid_ = self._parse_uuid(metadata['UUID'])
            provenance = self._parse_provenance(metadata['provenance'])
            type_ = self._parse_type(metadata['imports'], metadata['type'])
        return type_, uuid_, provenance

    # TODO this is mostly duplicated from WorkflowTemplate._parse_type.
    # Refactor!
    def _parse_type(self, imports, type_exp):
        type_exp = type_exp.split('\n')
        if len(type_exp) != 1:
            raise TypeError("Multiple lines in type expression of"
                            " artifact. Will not load to avoid arbitrary"
                            " code execution.")
        type_exp, = type_exp

        if ';' in type_exp:
            raise TypeError("Invalid type expression in artifact. Will not"
                            " load to avoid arbitrary code execution.")

        locals_ = {}
        for import_ in imports:
            path, class_ = import_.split(":")
            try:
                module = importlib.import_module(path)
            except ImportError:
                raise ImportError("The plugin which defines: %r is not"
                                  " installed." % path)
            class_ = getattr(module, class_)
            if not issubclass(class_, Type):
                raise TypeError("Non-Type artifact. Will not load to avoid"
                                " arbitrary code execution.")
            if class_.__name__ in locals_:
                raise TypeError("Duplicate type name (%r) in expression."
                                % class_.__name__)
            locals_[class_.__name__] = class_
        type_ = eval(type_exp, {'__builtins__': {}}, locals_)
        if not type_().is_concrete():
            raise TypeError("%r is not a concrete type. Only concrete types "
                            "can be loaded." % type_)
        return type_

    def _parse_uuid(self, string):
        return uuid.UUID(hex=string)

    def _parse_provenance(self, string):
        return string


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
