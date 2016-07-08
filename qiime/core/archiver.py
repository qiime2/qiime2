# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import os
import os.path
import shutil
import tempfile
import uuid
import yaml
import zipfile

import qiime.sdk
from .util import parse_type


# TODO use utf-8 encoding when reading/writing files
class Archiver:
    # This class will likely defer to one of many ArchiveFormats in the future.
    # There's only one supported archive format currently.
    _VERSION = '0.1.0'
    _VERSION_FILENAME = 'VERSION'
    _README_FILENAME = 'README.md'
    _METADATA_FILENAME = 'metadata.yaml'
    DATA_DIRNAME = 'data'

    @classmethod
    def extract(cls, filepath, output_dir):
        with zipfile.ZipFile(filepath, mode='r') as zf:
            zf.extractall(output_dir)
        return output_dir

    @classmethod
    def peek(cls, filepath):
        """

        Returns
        -------
        tuple
            Tuple of UUID, type, and provenance.

        """
        if not zipfile.is_zipfile(filepath):
            raise zipfile.BadZipFile(
                "%r is not a readable ZIP file, or the file does not exist" %
                filepath)

        root_dir = cls._get_root_dir(filepath)
        with zipfile.ZipFile(filepath, mode='r') as zf:
            version = cls._load_version(zf, root_dir)
            if version != cls._VERSION:
                raise ValueError(
                    "Unsupported archive format version %r. "
                    "Supported version(s): %r" % (version, cls._VERSION))

            return cls._load_metadata(zf, root_dir)

    @classmethod
    def load(cls, filepath):
        metadata = cls.peek(filepath)
        return cls(*metadata, archive_filepath=filepath)

    @classmethod
    def _get_root_dir(cls, filepath):
        return os.path.splitext(os.path.basename(filepath))[0]

    @classmethod
    def _load_version(cls, zf, root_dir):
        version_path = os.path.join(root_dir, cls._VERSION_FILENAME)
        with zf.open(version_path) as bytes_fh:
            with io.TextIOWrapper(bytes_fh, newline=None,
                                  encoding='utf-8') as fh:
                return fh.read().rstrip('\n')

    @classmethod
    def _load_metadata(cls, zf, root_dir):
        metadata_path = os.path.join(root_dir, cls._METADATA_FILENAME)
        with zf.open(metadata_path) as bytes_fh:
            with io.TextIOWrapper(bytes_fh, newline=None,
                                  encoding='utf-8') as fh:
                metadata = yaml.safe_load(fh)

        uuid_ = cls._parse_uuid(metadata['uuid'])
        type_ = parse_type(metadata['type'])
        provenance = cls._parse_provenance(metadata['provenance'])
        return uuid_, type_, provenance

    @classmethod
    def _parse_uuid(cls, string):
        return uuid.UUID(hex=string)

    # TODO implement provenance parsing for real
    @classmethod
    def _parse_provenance(cls, provenance):
        if isinstance(provenance, str):
            return None
        else:
            return qiime.sdk.Provenance(
                execution_uuid=uuid.UUID(hex=provenance['execution-uuid']),
                executor_reference=provenance['executor-reference'],
                artifact_uuids={k: uuid.UUID(hex=v) for k, v in
                                provenance['artifact-uuids'].items()},
                parameter_references=provenance['parameter-references']
            )

    def __init__(self, uuid, type, provenance, archive_filepath=None,
                 data_initializer=None):
        """

        Parameters
        ----------
        data_initializer : callable, optional
            Callable accepting a single str argument specifying the data
            directory to write data to. Callable should not return anything,
            return values will be ignored.

        """
        if archive_filepath is None and data_initializer is None:
            raise ValueError(
                "`archive_filepath` or `data_initializer` must be provided.")

        if archive_filepath is not None and data_initializer is not None:
            raise ValueError(
                "`archive_filepath` and `data_initializer` cannot both be "
                "provided.")

        self._uuid = uuid
        self._type = type
        self._provenance = provenance

        self._temp_dir = tempfile.mkdtemp(prefix='qiime2-archive-temp-')
        self._data_dir = os.path.join(self._temp_dir, self.DATA_DIRNAME)
        self._pid = os.getpid()

        if archive_filepath is not None:
            self._extract_data(archive_filepath)
        else:
            os.mkdir(self._data_dir)
            data_initializer(self._data_dir)

    def __del__(self):
        # Destructor can be called more than once.
        if os.path.exists(self._temp_dir) and self._pid == os.getpid():
            shutil.rmtree(self._temp_dir)

    def _extract_data(self, archive_filepath):
        root_dir = self._get_root_dir(archive_filepath)
        prefix = os.path.join(root_dir, self.DATA_DIRNAME, '')

        with zipfile.ZipFile(archive_filepath, mode='r') as zf:
            for file_ in zf.namelist():
                if file_.startswith(prefix) and file_ != prefix:
                    zf.extract(file_, path=self._temp_dir)
            shutil.move(os.path.join(self._temp_dir, prefix), self._temp_dir)
            os.rmdir(os.path.join(self._temp_dir, root_dir))

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

    def get_data_paths(self, recursive=True):
        iterable = iter(os.walk(self._data_dir))
        if not recursive:
            iterable = [next(iterable)]
        for root, _, files in iterable:
            for fp in files:
                fullpath = os.path.join(root, fp)
                yield (os.path.relpath(fullpath, self._data_dir), fullpath)

    def load_data(self, loader):
        return loader(self._data_dir)

    def save(self, filepath):
        root_dir = self._get_root_dir(filepath)

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
                        self.DATA_DIRNAME,
                        relpath
                    )
                    zf.write(abspath, arcname=archive_path)

    def _save_version(self, zf, root_dir):
        zf.writestr(os.path.join(root_dir, self._VERSION_FILENAME),
                    '%s\n' % self._VERSION)

    def _save_readme(self, zf, root_dir):
        zf.writestr(os.path.join(root_dir, self._README_FILENAME),
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
        zf.writestr(os.path.join(root_dir, self._METADATA_FILENAME),
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
                'executor-reference': self.provenance.executor_reference,
                'artifact-uuids': {k: str(v) for k, v in
                                   self.provenance.artifact_uuids.items()},
                'parameter-references': self.provenance.parameter_references
            }


_README_TEXT = """# QIIME 2 archive

This archive stores the data and associated metadata of a serialized QIIME 2
artifact or visualization.

**WARNING:** This is a temporary file format used for prototyping QIIME 2. Do
not rely on this format for any reason as it will not be supported past the
prototyping stage in its current form.
"""
