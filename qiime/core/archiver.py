# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import datetime
import io
import os
import tempfile
import time
import uuid
import yaml
import zipfile

import qiime.sdk


# TODO use utf-8 encoding when reading/writing files
class Archiver:
    # This class will likely defer to one of many ArchiveFormats in the future.
    # There's only one supported archive format currently.
    _version = '0.1.0'
    _version_path = 'VERSION'
    _readme_path = 'README.md'
    _metadata_path = 'metadata.yaml'
    _data_dir = 'data'

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
                    "Unsupported archive format version %r. "
                    "Supported version(s): %r" % (version, cls._version))

            uuid_, type_, provenance = cls._load_metadata(zf, root_dir)

        return cls(uuid_, type_, provenance, archive_filepath=filepath)

    @classmethod
    def _get_root_dir(cls, filepath):
        return os.path.splitext(os.path.basename(filepath))[0]

    @classmethod
    def _create_temp_dir(cls):
        return tempfile.TemporaryDirectory(prefix='qiime2-archive-temp-data-')

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
        type_ = cls._parse_type(metadata['type'])
        provenance = cls._parse_provenance(metadata['provenance'])
        return uuid_, type_, provenance

    @classmethod
    def _parse_uuid(cls, string):
        return uuid.UUID(hex=string)

    @classmethod
    def _parse_type(cls, type_exp):
        # TODO this is mostly duplicated from Workflow._parse_semantic_type.
        type_exp = type_exp.split('\n')
        if len(type_exp) != 1:
            raise TypeError(
                "Found multiple lines in archive type expression. Will not "
                "evaluate to avoid arbitrary code execution.")
        type_exp, = type_exp

        if ';' in type_exp:
            raise TypeError(
                "Invalid archive type expression. Will not evaluate to avoid "
                "arbitrary code execution.")

        pm = qiime.sdk.PluginManager()
        locals_ = {k: v[1] for k, v in pm.semantic_types.items()}
        locals_[qiime.core.type.Visualization.name] = \
            qiime.core.type.Visualization

        # TODO better error handling for un-eval-able type expressions
        return eval(type_exp, {'__builtins__': {}}, locals_)

    # TODO implement provenance parsing for real
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

    def __init__(self, uuid_, type_, provenance, archive_filepath=None):
        self._uuid = uuid_

        if isinstance(type_, str):
            type_ = self._parse_type(type_)
        self._type = type_

        self._provenance = provenance
        self._archive_filepath = archive_filepath
        self._temp_dir = None

    def __del__(self):
        if self._temp_dir is not None:
            self._temp_dir.cleanup()

    @property
    def uuid(self):
        return self._uuid

    @property
    def type(self):
        return self._type

    @property
    def provenance(self):
        return self._provenance

    def load_data(self, loader):
        if self._temp_dir is None:
            self._extract_data()

        return loader(self._temp_dir.name)

    def _extract_data(self):
        assert self._temp_dir is None and self._archive_filepath is not None

        self._temp_dir = self._create_temp_dir()
        temp_dir = self._temp_dir.name

        # TODO is there a better way to extract the files in an archive
        # subdirectory into a directory without keeping the nested archive
        # directory structure?
        with zipfile.ZipFile(self._archive_filepath, mode='r') as zf:
            root_dir = self._get_root_dir(self._archive_filepath)
            prefix = os.path.join(root_dir, self._data_dir, '')

            for file_ in zf.namelist():
                if file_.startswith(prefix):
                    extract_path = zf.extract(file_, path=temp_dir)
                    os.rename(extract_path,
                              os.path.join(temp_dir, os.path.basename(file_)))
            os.rmdir(os.path.join(temp_dir, prefix))

    def save_data(self, data, saver):
        assert self._temp_dir is None and self._archive_filepath is None

        self._temp_dir = self._create_temp_dir()
        saver(data, self._temp_dir.name)

    def save(self, filepath):
        root_dir = self._get_root_dir(filepath)

        if self._temp_dir is None:
            self._extract_data()

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
                        self._data_dir,
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
            'type': repr(self.type),
            'provenance': self._formatted_provenance()
        })
        zf.writestr(os.path.join(root_dir, self._metadata_path),
                    metadata_bytes)

    def _formatted_uuid(self):
        return str(self.uuid)

    def _formatted_provenance(self):
        if self.provenance is None:
            return ('No provenance available because this archive was '
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


_README_TEXT = """# QIIME 2 archive

This archive stores the data and associated metadata of a serialized QIIME 2
artifact or visualization.

**WARNING:** This is a temporary file format used for prototyping QIIME 2. Do
not rely on this format for any reason as it will not be supported past the
prototyping stage in its current form.
"""
