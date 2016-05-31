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
import tempfile
import yaml
import time
import datetime
import zipfile

import qiime.plugin
import qiime.sdk


# TODO use utf-8 encoding when reading/writing files
class Artifact:
    _version = '0.1.0'  # artifact archive format version
    _version_path = 'VERSION'
    _readme_path = 'README.md'
    _metadata_path = 'metadata.yaml'
    _archive_dir = 'data'

    @property
    def type(self):
        return self._type

    @property
    def provenance(self):
        return self._provenance

    @property
    def uuid(self):
        return self._uuid

    @classmethod
    def _get_root_dir(cls, filepath):
        return os.path.splitext(os.path.basename(filepath))[0]

    @classmethod
    def _create_temp_dir(cls):
        return tempfile.TemporaryDirectory(prefix='qiime2-temp-artifact-data-')

    @classmethod
    def load(cls, filepath):
        if not zipfile.is_zipfile(filepath):
            raise zipfile.BadZipFile(
                "%r is not a readable ZIP file, or the file does not exist" %
                filepath)

        root_dir = cls._get_root_dir(filepath)
        with zipfile.ZipFile(filepath, mode='r') as zf:
            version = cls._load_version(zf, root_dir)
            if version != cls._version:
                raise ValueError(
                    "Unsupported artifact archive format version %r. "
                    "Supported version(s): %r" % (version, cls._version))

            type_, uuid_, provenance = cls._load_metadata(zf, root_dir)

        artifact = cls.__new__(cls)

        artifact._type = type_
        artifact._provenance = provenance
        artifact._uuid = uuid_
        artifact._zip_filepath = filepath

        return artifact

    @classmethod
    def _load_version(cls, zf, root_dir):
        with zf.open(os.path.join(root_dir, cls._version_path)) as bytes_fh:
            with io.TextIOWrapper(bytes_fh, newline=None,
                                  encoding='utf-8') as fh:
                return fh.read().rstrip('\n')

    @classmethod
    def _load_metadata(cls, zf, root_dir):
        with zf.open(os.path.join(root_dir, cls._metadata_path)) as bytes_fh:
            with io.TextIOWrapper(bytes_fh, newline=None,
                                  encoding='utf-8') as fh:
                metadata = yaml.safe_load(fh)

        uuid_ = cls._parse_uuid(metadata['UUID'])
        provenance = cls._parse_provenance(metadata['provenance'])
        type_ = cls._parse_semantic_type(metadata['type'])
        return type_, uuid_, provenance

    # TODO this is mostly duplicated from Workflow._parse_semantic_type.
    @classmethod
    def _parse_semantic_type(cls, type_exp):
        type_exp = type_exp.split('\n')
        if len(type_exp) != 1:
            raise TypeError("Multiple lines in type expression of"
                            " artifact. Will not load to avoid arbitrary"
                            " code execution.")
        type_exp, = type_exp

        if ';' in type_exp:
            raise TypeError("Invalid type expression in artifact. Will not"
                            " load to avoid arbitrary code execution.")

        pm = qiime.sdk.PluginManager()
        locals_ = {k: v[1] for k, v in pm.semantic_types.items()}

        type_ = eval(type_exp, {'__builtins__': {}}, locals_)
        if not type_.is_concrete():
            raise TypeError("%r is not a concrete type. Only concrete types "
                            "can be loaded." % type_)
        return type_

    @classmethod
    def _parse_uuid(cls, string):
        return uuid.UUID(hex=string)

    # TODO implement provenance parsing
    @classmethod
    def _parse_provenance(cls, provenance):
        if isinstance(provenance, str):
            return None
        else:
            return qiime.sdk.Provenance(
                job_uuid=provenance['job-UUID'],
                artifact_uuids=provenance['artifact-uuids'],
                parameters=provenance['parameters'],
                workflow_reference=provenance['workflow-reference']
            )

    @classmethod
    def _from_view(cls, view, type_, provenance):
        """

        Parameters
        ----------
        view : Python object
            View to serialize.
        type_ : qiime.plugin.Type or str
            Semantic type of the artifact.
        provenance : qiime.sdk.Provenance
            Artifact provenance.

        """
        # TODO do some validation here, such as making sure `view` can be
        # written to `type_` archive format.
        if isinstance(type_, str):
            type_ = cls._parse_semantic_type(type_)
        if not type_.is_concrete():
            raise TypeError(
                "%r is not a concrete type. Artifacts can only contain "
                "concrete types." % type_)

        artifact = cls.__new__(cls)

        artifact._type = type_
        artifact._provenance = provenance
        artifact._uuid = uuid.uuid4()
        artifact._temp_dir = cls._create_temp_dir()
        artifact._save_view(view)

        return artifact

    def _save_view(self, view):
        # TODO update when types are ready
        type_format_map = {
            'FeatureTable[Frequency]': ('feature-table', 1),
            'FeatureTable[RelativeFrequency]': ('feature-table', 1),
            'FeatureTable[PresenceAbsence]': ('feature-table', 1),
            'TestType': ('test-archive-format', 1),
        }

        archive_format = type_format_map[repr(self.type)]
        archive_format = qiime.plugin.get_archive_format(*archive_format)

        view_type = type(view)
        writer = archive_format.writers[view_type]
        writer(view, self._temp_dir.name)

    def __init__(self):
        """

        This constructor is private.

        """
        raise NotImplementedError(
            "Artifact constructor is private. Use `Artifact.load` to "
            "construct an Artifact.")

    def __new__(cls):
        artifact = object.__new__(cls)

        artifact._type = None
        artifact._provenance = None
        artifact._uuid = None
        artifact._temp_dir = None
        artifact._zip_filepath = None

        return artifact

    def __del__(self):
        if self._temp_dir is not None:
            self._temp_dir.cleanup()

    def view(self, view_type):
        # TODO update when types are ready
        type_format_map = {
            'FeatureTable[Frequency]': ('feature-table', 1),
            'FeatureTable[RelativeFrequency]': ('feature-table', 1),
            'FeatureTable[PresenceAbsence]': ('feature-table', 1),
            'TestType': ('test-archive-format', 1),
        }

        archive_format = type_format_map[repr(self.type)]
        archive_format = qiime.plugin.get_archive_format(*archive_format)

        reader = archive_format.readers[view_type]

        if self._temp_dir is None:
            self._extract_artifact_data()

        return reader(self._temp_dir.name)

    def _extract_artifact_data(self):
        assert self._temp_dir is None and self._zip_filepath is not None

        self._temp_dir = self._create_temp_dir()
        temp_dir = self._temp_dir.name

        # TODO is there a better way to extract the files in an archive
        # subdirectory into a directory without keeping the nested archive
        # directory structure?
        with zipfile.ZipFile(self._zip_filepath, mode='r') as zf:
            root_dir = self._get_root_dir(self._zip_filepath)
            prefix = os.path.join(root_dir, self._archive_dir, '')

            for file_ in zf.namelist():
                if file_.startswith(prefix):
                    extract_path = zf.extract(file_, path=temp_dir)
                    os.rename(extract_path,
                              os.path.join(temp_dir, os.path.basename(file_)))
            os.rmdir(os.path.join(temp_dir, prefix))

    def save(self, filepath):
        root_dir = self._get_root_dir(filepath)

        if self._temp_dir is None:
            self._extract_artifact_data()

        with zipfile.ZipFile(filepath, mode='w',
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            self._save_version(zf, root_dir)
            self._save_readme(zf, root_dir)
            self._save_metadata(zf, root_dir)

            for root, dirs, files in os.walk(self._temp_dir.name):
                for file_ in files:
                    temp_path = os.path.join(root, file_)
                    archive_path = os.path.join(
                        root_dir,
                        self._archive_dir,
                        root.replace(self._temp_dir.name, '', 1),
                        file_
                    )
                    zf.write(temp_path, arcname=archive_path)

    def _save_version(self, zf, root_dir):
        zf.writestr(os.path.join(root_dir, self._version_path),
                    '%s\n' % self._version)

    def _save_readme(self, zf, root_dir):
        zf.writestr(os.path.join(root_dir, self._readme_path),
                    _README_TEXT)

    # TODO clean up metadata yaml formatting. It currently dumps Python
    # objects, `yaml.safe_dump` call needs to be updated to format lists and
    # dicts as typical yaml.
    def _save_metadata(self, zf, root_dir):
        metadata_bytes = yaml.safe_dump({
            'UUID': self._formatted_uuid(),
            'provenance': self._formatted_provenance(),
            'type': repr(self.type)
        })
        zf.writestr(os.path.join(root_dir, self._metadata_path),
                    metadata_bytes)

    def _formatted_uuid(self):
        return str(self.uuid)

    def _formatted_provenance(self):
        if self.provenance is None:
            return ('No provenance available because this artifact was '
                    'generated independently of QIIME.')
        else:
            return {
                'job-UUID': str(self.provenance.job_uuid),
                'artifact-uuids': self.provenance.artifact_uuids,
                # TODO: need to store parameter types (or is it enough
                # that the workflow template reference could get us that?)
                'parameters': self.provenance.parameters,
                'runtime': self._formatted_runtime(),
                'workflow-reference': self.provenance.workflow_reference
            }

    def _formatted_runtime(self):
        stop_time = time.time()
        run_time = stop_time - self.provenance.start_time
        datetime_start = datetime.datetime.fromtimestamp(
            self.provenance.start_time)
        datetime_stop = datetime.datetime.fromtimestamp(stop_time)
        return {
            'run-time-seconds': run_time,
            'execution-start-time': datetime_start.isoformat(),
            'execution-stop-time': datetime_stop.isoformat()
        }


_README_TEXT = """# QIIME 2 artifact

This archive stores the data and associated metadata of a serialized QIIME 2
artifact.

**WARNING:** This is a temporary file format used for prototyping QIIME 2. Do
not rely on this format for any reason as it will not be supported past the
prototyping stage in its current form.
"""
