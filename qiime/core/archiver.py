# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import os
import shutil
import tempfile
import uuid
import yaml
import zipfile

import qiime.sdk


# TODO use utf-8 encoding when reading/writing files
class Archiver:
    # This class will likely defer to one of many ArchiveFormats in the future.
    # There's only one supported archive format currently.
    _version = '0.1.0'
    _version_filename = 'VERSION'
    _readme_filename = 'README.md'
    _metadata_filename = 'metadata.yaml'
    _data_dirname = 'data'

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
    def _load_version(cls, zf, root_dir):
        version_path = os.path.join(root_dir, cls._version_filename)
        with zf.open(version_path) as bytes_fh:
            with io.TextIOWrapper(bytes_fh, newline=None,
                                  encoding='utf-8') as fh:
                return fh.read().rstrip('\n')

    @classmethod
    def _load_metadata(cls, zf, root_dir):
        metadata_path = os.path.join(root_dir, cls._metadata_filename)
        with zf.open(metadata_path) as bytes_fh:
            with io.TextIOWrapper(bytes_fh, newline=None,
                                  encoding='utf-8') as fh:
                metadata = yaml.safe_load(fh)

        uuid_ = cls._parse_uuid(metadata['uuid'])
        type_ = cls._parse_type(metadata['type'])
        provenance = cls._parse_provenance(metadata['provenance'])
        return uuid_, type_, provenance

    @classmethod
    def _parse_uuid(cls, string):
        return uuid.UUID(hex=string)

    @classmethod
    def _parse_type(cls, type_exp):
        # TODO this is mostly duplicated from Method._parse_semantic_type,
        # refactor.
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
        locals_['Properties'] = qiime.core.type.semantic.Properties

        # TODO better error handling for un-eval-able type expressions
        return eval(type_exp, {'__builtins__': {}}, locals_)

    # TODO implement provenance parsing for real
    @classmethod
    def _parse_provenance(cls, provenance):
        if isinstance(provenance, str):
            return None
        else:
            return qiime.sdk.Provenance(
                execution_uuid=provenance['execution-uuid'],
                method_reference=provenance['method-reference'],
                artifact_uuids=provenance['artifact-uuids'],
                parameters=provenance['parameters']
            )

    def __init__(self, uuid, type, provenance, archive_filepath=None):
        self._uuid = uuid

        if isinstance(type, str):
            type = self._parse_type(type)
        self._type = type

        self._provenance = provenance
        self._archive_filepath = archive_filepath

        self._temp_dir = tempfile.mkdtemp(prefix='qiime2-archive-temp-')
        self._data_dir = os.path.join(self._temp_dir, self._data_dirname)
        self._pid = os.getpid()

    def __del__(self):
        # Destructor can be called more than once.
        if os.path.exists(self._temp_dir) and self._pid == os.getpid():
            shutil.rmtree(self._temp_dir)

    @property
    def uuid(self):
        return self._uuid

    @property
    def type(self):
        return self._type

    @property
    def provenance(self):
        return self._provenance

    def orphan(self, pid):
        self._pid = pid

    def load_data(self, loader):
        if not os.path.exists(self._data_dir):
            self._extract_data()

        return loader(self._data_dir)

    def _extract_data(self):
        assert not os.path.exists(self._data_dir)
        assert self._archive_filepath is not None

        with zipfile.ZipFile(self._archive_filepath, mode='r') as zf:
            root_dir = self._get_root_dir(self._archive_filepath)
            prefix = os.path.join(root_dir, self._data_dirname, '')

            for file_ in zf.namelist():
                if file_.startswith(prefix) and file_ != prefix:
                    zf.extract(file_, path=self._temp_dir)
            shutil.move(os.path.join(self._temp_dir, prefix), self._temp_dir)
            os.rmdir(os.path.join(self._temp_dir, root_dir))

    def save_data(self, data, saver):
        assert not os.path.exists(self._data_dir)
        assert self._archive_filepath is None

        os.mkdir(self._data_dir)
        saver(data, self._data_dir)

    def save(self, filepath):
        root_dir = self._get_root_dir(filepath)

        if not os.path.exists(self._data_dir):
            self._extract_data()

        with zipfile.ZipFile(filepath, mode='w',
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            self._save_version(zf, root_dir)
            self._save_readme(zf, root_dir)
            self._save_metadata(zf, root_dir)

            for root, dirs, files in os.walk(self._data_dir):
                for file_ in files:
                    abspath = os.path.join(root, file_)
                    relpath = os.path.relpath(abspath, start=self._data_dir)
                    archive_path = os.path.join(
                        root_dir,
                        self._data_dirname,
                        relpath
                    )
                    zf.write(abspath, arcname=archive_path)

    def _save_version(self, zf, root_dir):
        zf.writestr(os.path.join(root_dir, self._version_filename),
                    '%s\n' % self._version)

    def _save_readme(self, zf, root_dir):
        zf.writestr(os.path.join(root_dir, self._readme_filename),
                    _README_TEXT)

    # TODO clean up metadata yaml formatting. It currently dumps Python
    # objects, `yaml.safe_dump` call needs to be updated to format lists and
    # dicts as typical yaml.
    def _save_metadata(self, zf, root_dir):
        metadata_bytes = yaml.safe_dump({
            'uuid': self._formatted_uuid(),
            'type': repr(self.type),
            'provenance': self._formatted_provenance()
        })
        zf.writestr(os.path.join(root_dir, self._metadata_filename),
                    metadata_bytes)

    def _formatted_uuid(self):
        return str(self.uuid)

    # TODO this is a provenance placeholder for now
    def _formatted_provenance(self):
        if self.provenance is None:
            return ('No provenance available because this archive was '
                    'generated independently of QIIME.')
        else:
            return {
                'execution-uuid': str(self.provenance.execution_uuid),
                'method-reference': self.provenance.method_reference,
                'artifact-uuids': self.provenance.artifact_uuids,
                'parameters': {
                    k: str(v) for k, v in self.provenance.parameters.items()}
            }


_README_TEXT = """# QIIME 2 archive

This archive stores the data and associated metadata of a serialized QIIME 2
artifact or visualization.

**WARNING:** This is a temporary file format used for prototyping QIIME 2. Do
not rely on this format for any reason as it will not be supported past the
prototyping stage in its current form.
"""
