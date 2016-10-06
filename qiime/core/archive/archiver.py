import collections
import uuid as _uuid
import pathlib
import zipfile
import tempfile
import importlib
import os
import shutil
import io

import qiime

_VERSION_TEMPLATE = """\
QIIME 2
archive: %s
framework: %s
"""

ArchiveRecord = collections.namedtuple(
    'ArchiveRecord', ['root', 'version_fp', 'uuid', 'version',
                      'framework_version'])


class _Archive:
    """Abstraction layer over the archive filesystem.

    Responsible for details concerning manipulating an archive agnostic to its
    format. It is responsible for managing archive UUID, format, and framework
    versions as those are designed to be constant throughout all future
    format implementations. Breaking compatibility with that is a BIG DEAL and
    should avoided at (nearly) any cost.

    Example filesystem::

        <archive root>/
        !--- 770509e6-85f4-432c-9663-cdc04eb07db2
            |--- VERSION
            !--- <whatever format defines>

    VERSION file::

        QIIME 2
        archive: <archive version>
        framework: <framework version>

    This file is itentionally not YAML/INI/An actual format. This is to
    discourage the situation where the format changes from something like YAML
    to another format and VERSION is updated with it "for consistency".

    To emphasize, the VERSION (filepath and content) and root archive structure
    MUST NOT CHANGE. If they change, then there is no longer a consistent way
    to dispatch to an appropriate format.

    """
    VERSION_FILE = 'VERSION'

    @classmethod
    def is_archive_type(cls, filepath):
        raise NotImplementedError

    @classmethod
    def setup(cls, path, version, framework_version):
        uuid = _uuid.uuid4()


        root_dir = path / str(uuid)
        root_dir.mkdir()
        version_fp = root_dir / cls.VERSION_FILE

        version_fp.write_text(_VERSION_TEMPLATE % (version, framework_version))

        return ArchiveRecord(root_dir, version_fp, uuid, version,
                             framework_version)

    @classmethod
    def save(cls, source, destination):
        raise NotImplementedError

    def __init__(self, path):
        self.path = path

        self.uuid = self._get_uuid()
        self.version, self.framework_version = self._get_versions()

    def _get_uuid(self):
        if not self.path.exists():
            raise TypeError("%s does not exist or is not a filepath."
                            % self.path)

        roots = set()
        for relpath in self.relative_iterdir():
            if not relpath.startswith('.'):
                roots.add(relpath)

        if len(roots) == 0:
            raise ValueError("Archive does not have a visible root directory.")
        if len(roots) > 1:
            raise ValueError("Archive has multiple root directories: %r"
                             % roots)
        uuid = roots.pop()
        if not self._is_uuid4(uuid):
            raise ValueError(
                "Archive root directory name %r is not a valid version 4 "
                "UUID." % root)
        return uuid

    def _is_uuid4(self, uuid_str):
        # Adapted from https://gist.github.com/ShawnMilo/7777304
        try:
            uuid = _uuid.UUID(hex=uuid_str, version=4)
        except ValueError:
            # The string is not a valid hex code for a UUID.
            return False

        # If uuid_str is a valid hex code, but an invalid uuid4, UUID.__init__
        # will convert it to a valid uuid4.
        return str(uuid) == uuid_str

    def _get_versions(self):
        try:
            with self.open(self.VERSION_FILE) as fh:
                header, version_line, framework_version_line, eof = \
                    fh.read().split('\n')
            assert(header.strip() == 'QIIME 2')
            version = version_line.split(':')[1].strip()
            framework_version = framework_version_line.split(':')[1].strip()
            return version, framework_version
        except Exception:
            # TODO: make a "better" parser which isn't just a catch-all
            raise ValueError("Archive does not contain a correctly formatted"
                             " VERSION file.")

    def relative_iterdir(self, relpath='.'):
        raise NotImplementedError

    def open(self, relpath):
        raise NotImplementedError

    def mount(self, filepath):
        raise NotImplementedError



class _ZipArchive(_Archive):
    """A specific variant of Archive which deals with ZIP64 files."""

    @classmethod
    def is_archive_type(cls, path):
        return zipfile.is_zipfile(str(path))

    @classmethod
    def extract(self, source, destination):
        with zipfile.ZipFile(str(source), mode='r') as zf:
            zf.extractall(str(destination))

    @classmethod
    def save(cls, source, destination):
        with zipfile.ZipFile(destination, mode='w',
                             compression=zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zf:
            for root, dirs, files in os.walk(str(source)):
                # Prune hidden directories from traversal. Strategy modified
                # from http://stackoverflow.com/a/13454267/3776794
                dirs[:] = [d for d in dirs if not d.startswith('.')]

                for file in files:
                    if file.startswith('.'):
                        continue

                    abspath = pathlib.Path(root) / file
                    relpath = abspath.relative_to(source)

                    zf.write(str(abspath), arcname=cls._as_zip_path(relpath))

    def relative_iterdir(self, relpath='.'):
        relpath = self._as_zip_path(relpath)
        relpath_len = len(relpath)
        seen = set()
        with zipfile.ZipFile(str(self.path), mode='r') as zf:
            for name in zf.namelist():
                if name.startswith(relpath):
                    parts = pathlib.PurePosixPath(name).parts
                    if len(parts) > 0:
                        result = parts[0]
                        if result not in seen:
                            seen.add(result)
                            yield result

    def open(self, relpath):
        relpath = pathlib.Path(str(self.uuid)) / relpath
        with zipfile.ZipFile(str(self.path), mode='r') as zf:
            # The filehandle will still work even when `zf` is "closed"
            return io.TextIOWrapper(zf.open(self._as_zip_path(relpath)))

    def mount(self, filepath):
        # TODO: use FUSE/MacFUSE/Dokany bindings (many Python bindings are
        # outdated, we may need to take up maintenance/fork)
        root = self.extract(filepath)
        return ArchiveRecord(root, root / self.VERSION_FILE,
                             self.uuid, self.version, self.framework_version)

    def extract(self, filepath):
        with zipfile.ZipFile(str(self.path), mode='r') as zf:
            zf.extractall(str(filepath))
        return pathlib.Path(filepath) / str(self.uuid)

    @classmethod
    def _as_zip_path(self, path):
        return str(pathlib.PurePosixPath(path))


class Archiver:
    CURRENT_FORMAT_VERSION = '1'
    CURRENT_ARCHIVE = _ZipArchive
    _FORMAT_REGISTRY = {
        # NOTE: add more archive formats as things change
        '1': 'qiime.core.archive.format.v1:ArchiveFormat'
    }

    @classmethod
    def _make_temp_path(cls):
        # Not using OutPath because it cleans itself up. Archiver should be
        # responsible for that.
        # TODO: normalize `mkdtemp` when we have framework temp locations.
        return pathlib.Path(tempfile.mkdtemp(prefix='qiime2-archive-'))

    @classmethod
    def get_format_class(cls, version):
        try:
            imp, fmt_cls = cls._FORMAT_REGISTRY[version].split(':')
        except KeyError:
            return None
        return getattr(importlib.import_module(imp), fmt_cls)

    @classmethod
    def get_archive(cls, filepath):
        filepath = pathlib.Path(filepath)
        if not filepath.exists():
            raise ValueError("%s does not exist." % filepath)

        if _ZipArchive.is_archive_type(filepath):
            archive = _ZipArchive(filepath)
        else:
            raise ValueError("%s is not a QIIME archive." % filepath)

        return archive

    @classmethod
    def _futuristic_archive_error(filepath, archive):
        raise ValueError("%s was created by 'QIIME %s'. The currently"
                         " installed framework cannot interpret archive"
                         " version %r."
                         % (filepath, archive.framework_version,
                            archive.version))

    @classmethod
    def peek(cls, filepath):
        archive = cls.get_archive(filepath)
        Format = cls.get_format_class(archive.version)
        if Format is None:
            cls._futuristic_archive_error(filepath, archive)
        # NOTE: in the future, we may want to manipulate the results so that
        # older formats provide the "new" API even if they don't support it.
        # e.g. a new format has a new property that peek should describe. We
        # add some compatability code here to return a default for that
        # property on older formats.
        return Format.load_metadata(archive)

    @classmethod
    def extract(cls, filepath, dest):
        archive = cls.get_archive(filepath)
        # Format really doesn't matter, the archive knows how to extract so
        # that is sufficient, furthermore it would suck if something was wrong
        # with an archive's format and extract failed to actually extract.
        return archive.extract(dest)

    @classmethod
    def load(cls, filepath):
        archive = cls.get_archive(filepath)
        Format = cls.get_format_class(archive.version)
        if Format is None:
            cls._futuristic_archive_error(filepath, archive)

        path = cls._make_temp_path()
        rec = archive.mount(path)

        return cls(path, Format(rec))

    @classmethod
    def from_data(cls, type, format, data_initializer, provenance_capture):
        path = cls._make_temp_path()
        rec = cls.CURRENT_ARCHIVE.setup(path, cls.CURRENT_FORMAT_VERSION,
                                        qiime.__version__)

        Format = cls.get_format_class(cls.CURRENT_FORMAT_VERSION)
        Format.write(rec, type, format, data_initializer, provenance_capture)

        return cls(path, Format(rec))

    def __init__(self, path, fmt):
        self.path = path
        self._fmt = fmt
        self._pid = os.getpid()

    @property
    def uuid(self):
        return self._fmt.uuid

    @property
    def type(self):
        return self._fmt.type

    @property
    def format(self):
        return self._fmt.format

    @property
    def provenance_dir(self):
        return self._fmt.provenance_dir

    @property
    def data_dir(self):
        return self._fmt.data_dir

    def __del__(self):
        # Destructor can be called more than once.
        if self.path.exists() and self._pid == os.getpid():
            shutil.rmtree(str(self.path))

    def save(self, filepath):
        self.CURRENT_ARCHIVE.save(self.path, filepath)

    def orphan(self, pid):
        self._pid = pid
