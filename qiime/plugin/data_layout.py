# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import os.path

from .file_format import FileFormat


class ValidationError(Exception):
    pass


class DataLayout:
    def __init__(self, name, version):
        self.name = name
        self.version = version
        # Relpath (str) mapped to FileFormat class. The code in this class
        # assumes that each relpath key represents a single file. glob/regular
        # expressions are not supported yet.
        self.files = {}
        self._finalized = False
        self._reader_views = {}
        self._writer_views = {}

    def register_file(self, relpath, file_format_type):
        if self._finalized:
            raise RuntimeError(
                "%r has been finalized and cannot have additional files "
                "registered to it." % self)
        if relpath in self.files:
            raise ValueError(
                "File %r has already been registered with this data layout."
                % relpath)
        if not issubclass(file_format_type, FileFormat):
            raise TypeError(
                "`file_format_type` must be a subclass of "
                "qiime.plugin.FileFormat.")
        self.files[relpath] = file_format_type

    def finalize(self):
        """Prevent data layout from having additional files registered on it.

        Call this method after all files defining this data layout have been
        registered with `register_file`.

        Data layout readers and writers can continue to be registered after
        finalization as those operations don't alter the definition of the data
        layout itself.

        This method is idempotent.

        """
        self._finalized = True

    def validate(self, path):
        if not os.path.exists(path):
            raise OSError(
                "Path does not exist or is not accessible: %r" % path)

        if os.path.isfile(path):
            if len(self.files) != 1:
                raise NotADirectoryError(
                    "%r requires %d files but only a single file %r was "
                    "provided. A directory must be provided to this data "
                    "layout." % (self, len(self.files), path))
            file_format = next(iter(self.files.values()))
            path_formats = [(path, file_format)]
        else:
            # `path` is a directory
            filepaths = []
            for root, dirs, files in os.walk(path):
                for dir in dirs:
                    dirpath = os.path.join(root, dir)
                    if len(os.listdir(dirpath)) == 0:
                        raise ValidationError(
                            "Found empty directory %r. Data layouts do "
                            "not support empty directories." % dirpath)
                for file in files:
                    filepaths.append(os.path.join(root, file))

            if len(filepaths) != len(self.files):
                raise ValidationError(
                    "%r requires %d files. Found %d files in %r. This data "
                    "layout requires the following filepaths: %s"
                    % (self, len(self.files), len(filepaths), path,
                       self._formatted_relpaths()))

            path_formats = []
            for filepath in filepaths:
                relpath = os.path.relpath(filepath, start=path)
                if relpath not in self.files:
                    raise ValidationError(
                        "File %r is not a filepath recognized by %r. This "
                        "data layout requires the following filepaths: %s"
                        % (filepath, self, self._formatted_relpaths()))
                file_format = self.files[relpath]
                path_formats.append((filepath, file_format))

        for filepath, file_format in path_formats:
            if not file_format.sniff(filepath):
                raise ValidationError(
                    "File %r does not appear to be a %r file."
                    % (filepath, file_format.name))

    def _formatted_relpaths(self):
        return ', '.join([str(relpath) for relpath in sorted(self.files)])

    def __repr__(self):
        return '%s(name=%r, version=%r)' % (self.__class__.__name__, self.name,
                                            self.version)

    @property
    def readers(self):
        return self._reader_views

    @property
    def writers(self):
        return self._writer_views

    def reader(self, view_type):
        def decorator(reader_function):
            self._reader_views[view_type] = reader_function
            return reader_function
        return decorator

    def writer(self, view_type):
        def decorator(writer_function):
            self._writer_views[view_type] = writer_function
            return writer_function
        return decorator
