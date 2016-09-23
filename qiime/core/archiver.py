# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import io
import os
import os.path
import shutil
import tempfile
import uuid
import yaml
import zipfile

import qiime.sdk

# Allow OrderedDict to be serialized for YAML representation
yaml.add_representer(collections.OrderedDict, lambda dumper, data:
                     dumper.represent_dict(data.items()))


class Archiver:
    # This class will likely defer to one of many ArchiveFormats in the future.
    # There's only one supported archive format currently.
    _VERSION = '0.3.0'
    _VERSION_FILENAME = 'VERSION'
    _README_FILENAME = 'README.md'
    _METADATA_FILENAME = 'metadata.yaml'
    DATA_DIRNAME = 'data'
    # Zipfiles always use forward slash as the path separator.
    _ZIPSEP = '/'

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
            Tuple of UUID, type, format, and provenance.

        """
        if not zipfile.is_zipfile(filepath):
            raise zipfile.BadZipFile(
                "%r is not a readable ZIP file, or the file does not exist" %
                filepath)

        with zipfile.ZipFile(filepath, mode='r') as zf:
            root_dir = cls._get_root_dir(zf)
            version = cls._load_version(zf, root_dir)
            if version != cls._VERSION:
                raise ValueError(
                    "Unsupported archive format version %r. "
                    "Supported version(s): %r" % (version, cls._VERSION))

            return cls._load_metadata(zf, root_dir)

    @classmethod
    def load(cls, filepath):
        def extract_data(data_dir):
            with zipfile.ZipFile(filepath, mode='r') as zf:
                root_dir = cls._get_root_dir(zf)

                # This archive path is <archive uuid>/data/
                # The trailing slash is necessary because directory entries
                # in the zipfile have trailing slashes, and these directory
                # entries shouldn't be extracted (see filtering logic in loop
                # below).
                archive_data_dir = cls._ZIPSEP.join(
                    [root_dir, cls.DATA_DIRNAME, ''])

                for file_ in zf.namelist():
                    if file_.startswith(archive_data_dir) and \
                            file_ != archive_data_dir:
                        zf.extract(file_, path=data_dir)

            # This filesystem path is `data_dir`/<archive uuid>/data/
            # Archive path separators must be converted to OS-specific
            # separators because this is a filesystem path now instead of an
            # archive path.
            extracted_data_dir = os.path.join(
                data_dir, *archive_data_dir.split(cls._ZIPSEP))

            for name in os.listdir(extracted_data_dir):
                # Note: if the archive's data directory contains a top-level
                # file or directory with a name matching this archive's UUID,
                # `shutil.move` will raise an error stating the destination
                # path already exists. This will likely never happen in
                # practice; it is similar to a UUID collision. If it does, the
                # user will get a somewhat cryptic error message, report it,
                # and we can make this code more robust. Until then, it isn't
                # worth complicating this code any further.
                shutil.move(os.path.join(extracted_data_dir, name), data_dir)

            # Now that all extracted data has been moved to `data_dir`,
            # recursively remove empty directories starting at the leaf.
            os.removedirs(extracted_data_dir)

        uuid_, type, format, provenance = cls.peek(filepath)
        return cls(type, format, provenance, data_initializer=extract_data,
                   uuid_=uuid_)

    @classmethod
    def _get_root_dir(cls, zf):
        roots = set()
        for relpath in zf.namelist():
            if not relpath.startswith('.'):
                root = relpath.split(cls._ZIPSEP, maxsplit=1)[0]
                roots.add(root)
        if len(roots) == 0:
            raise ValueError("Archive does not have a visible root directory.")
        if len(roots) > 1:
            raise ValueError(
                "Archive has multiple root directories: %r" % roots)

        root = roots.pop()
        if not cls._is_uuid4(root):
            raise ValueError(
                "Archive root directory name %r is not a valid version 4 "
                "UUID." % root)
        return root

    @classmethod
    def _is_uuid4(cls, uuid_str):
        # Adapted from https://gist.github.com/ShawnMilo/7777304
        try:
            uuid_ = uuid.UUID(uuid_str, version=4)
        except ValueError:
            # The string is not a valid hex code for a UUID.
            return False

        # If uuid_str is a valid hex code, but an invalid uuid4, UUID.__init__
        # will convert it to a valid uuid4.
        return str(uuid_) == uuid_str

    @classmethod
    def _load_version(cls, zf, root_dir):
        version_path = root_dir + cls._ZIPSEP + cls._VERSION_FILENAME
        with zf.open(version_path) as bytes_fh:
            with io.TextIOWrapper(bytes_fh, newline=None,
                                  encoding='utf-8', errors='strict') as fh:
                return fh.read().rstrip('\n')

    @classmethod
    def _load_metadata(cls, zf, root_dir):
        metadata_path = root_dir + cls._ZIPSEP + cls._METADATA_FILENAME
        with zf.open(metadata_path) as bytes_fh:
            with io.TextIOWrapper(bytes_fh, newline=None,
                                  encoding='utf-8', errors='strict') as fh:
                metadata = yaml.safe_load(fh)

        uuid_ = cls._parse_uuid(metadata['uuid'])
        if root_dir != str(uuid_):
            raise ValueError(
                "Archive root directory must match UUID present in archive's "
                "metadata: %r != %r" % (root_dir, str(uuid_)))

        type_ = qiime.sdk.parse_type(metadata['type'])
        provenance = cls._parse_provenance(metadata['provenance'])
        format = metadata['format']
        return uuid_, type_, format, provenance

    @classmethod
    def _parse_uuid(cls, string):
        if not cls._is_uuid4(string):
            raise ValueError("%r is not a valid version 4 UUID." % string)
        return uuid.UUID(hex=string, version=4)

    # TODO implement provenance parsing for real
    @classmethod
    def _parse_provenance(cls, provenance):
        if isinstance(provenance, str):
            if provenance == 'None':
                return None
            else:
                return provenance
        else:
            return qiime.sdk.Provenance(
                execution_uuid=cls._parse_uuid(provenance['execution-uuid']),
                # TODO this should be 'action-reference' to match QIIME 2
                # nomenclature. Not worth updating now while provenance is
                # stubbed.
                executor_reference=provenance['executor-reference'],
                artifact_uuids={k: cls._parse_uuid(v) for k, v in
                                provenance['artifact-uuids'].items()},
                parameter_references=provenance['parameter-references']
            )

    def __init__(self, type, format: str, provenance, data_initializer,
                 uuid_=None):
        """

        Parameters
        ----------
        data_initializer : callable
            Callable accepting a single str argument specifying the data
            directory to write data to. Callable should not return anything;
            return values will be ignored.
        uuid_ : uuid.UUID, optional
            Version 4 UUID. If not provided, a random version 4 UUID will be
            generated.

        """
        if uuid_ is None:
            uuid_ = uuid.uuid4()
        else:
            if not isinstance(uuid_, uuid.UUID) or uuid_.version != 4:
                raise TypeError(
                    "`uuid_` must be a version 4 UUID, not %r." % uuid_)
        self._uuid = uuid_

        self._type = type
        self._provenance = provenance
        self._format = format

        # TODO: set _temp_dir to be read-only, don't forget that windows will
        # fail on shutil.rmtree, so add an onerror callback which chmods and
        # removes again
        self._temp_dir = tempfile.mkdtemp(prefix='qiime2-archive-temp-')
        self._data_dir = os.path.join(self._temp_dir, self.DATA_DIRNAME)
        os.mkdir(self._data_dir)
        data_initializer(self._data_dir)

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

    @property
    def format(self):
        return self._format

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
        if filepath == '':
            raise ValueError("Cannot save to an empty filepath.")
        if os.path.isdir(filepath):
            raise ValueError("Cannot save to a directory (%s)." % filepath)
        if filepath.endswith(os.sep):
            raise ValueError(
                "Cannot save to a filepath ending in a path separator (%s)." %
                filepath)

        root_dir = self._formatted_uuid()
        with zipfile.ZipFile(filepath, mode='w',
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            self._save_version(zf, root_dir)
            self._save_readme(zf, root_dir)
            self._save_metadata(zf, root_dir)

            for root, dirs, files in os.walk(self._data_dir):
                # Prune hidden directories from traversal. Strategy modified
                # from http://stackoverflow.com/a/13454267/3776794
                dirs[:] = [d for d in dirs if not d.startswith('.')]

                for file_ in files:
                    if file_.startswith('.'):
                        continue

                    abspath = os.path.join(root, file_)
                    relpath = os.path.relpath(abspath, start=self._data_dir)
                    archive_path = self._ZIPSEP.join([
                        root_dir,
                        self.DATA_DIRNAME,
                        relpath.replace(os.sep, self._ZIPSEP)
                    ])
                    zf.write(abspath, arcname=archive_path)

    def _save_version(self, zf, root_dir):
        version = '%s\n' % self._VERSION
        zf.writestr(root_dir + self._ZIPSEP + self._VERSION_FILENAME,
                    version.encode('utf-8'))

    def _save_readme(self, zf, root_dir):
        zf.writestr(root_dir + self._ZIPSEP + self._README_FILENAME,
                    _README_TEXT.encode('utf-8'))

    def _save_metadata(self, zf, root_dir):
        metadata_bytes = yaml.dump(collections.OrderedDict([
            ('uuid', self._formatted_uuid()),
            ('type', repr(self.type)),
            ('format', self.format),
            ('provenance', self._formatted_provenance())
        ]), default_flow_style=False)
        zf.writestr(root_dir + self._ZIPSEP + self._METADATA_FILENAME,
                    metadata_bytes.encode('utf-8'))

    def _formatted_uuid(self):
        return str(self.uuid)

    # TODO this is a provenance placeholder for now
    def _formatted_provenance(self):
        if self.provenance is None:
            return str(None)
        elif isinstance(self.provenance, str):
            return self.provenance
        else:
            return {
                'execution-uuid': str(self.provenance.execution_uuid),
                # TODO this should be 'action-reference' to match QIIME 2
                # nomenclature. Not worth updating now while provenance is
                # stubbed.
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
