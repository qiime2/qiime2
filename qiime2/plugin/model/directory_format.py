# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import pathlib

from qiime2.core import transform
from .base import FormatBase, ValidationError, _check_validation_level


class PathMakerDescriptor:
    def __init__(self, file):
        self.file = file

    def __get__(self, obj, cls=None):
        if obj is None:
            raise Exception()
        return getattr(obj, self.file.name).path_maker


class File:
    def __init__(self, pathspec, *, format=None):
        if format is None:
            raise TypeError("Must provide a format.")
        self.pathspec = pathspec
        self.format = format

    def __get__(self, obj, cls=None):
        if obj is None:
            return self
        return BoundFile(self.name, self.pathspec, self.format, obj)


class FileCollection(File):
    def __init__(self, pathspec, *, format=None):
        super().__init__(pathspec, format=format)
        self._path_maker = None

    def set_path_maker(self, function):
        self._path_maker = function
        return PathMakerDescriptor(self)

    def __get__(self, obj, cls=None):
        if obj is None:
            return self

        if self._path_maker is None:
            raise NotImplementedError()

        return BoundFileCollection(self.name, self.pathspec, self.format,
                                   obj, path_maker=self._path_maker)


class BoundFile:
    @property
    def mode(self):
        return self._directory_format._mode

    def __init__(self, name, pathspec, format, directory_format):
        self.name = name
        self.pathspec = pathspec
        self.format = format
        self._directory_format = directory_format
        self._path_maker = lambda s: pathspec

    def view(self, view_type):
        from_type = transform.ModelType.from_view_type(self.format)
        to_type = transform.ModelType.from_view_type(view_type)

        transformation = from_type.make_transformation(to_type)
        return transformation(self.path_maker())

    def write_data(self, view, view_type, **kwargs):
        # TODO: make `view_type` optional like in `Artifact.import_data`
        if self.mode != 'w':
            raise TypeError("Cannot use `set`/`add` when mode=%r" % self.mode)
        from_type = transform.ModelType.from_view_type(view_type)
        to_type = transform.ModelType.from_view_type(self.format)

        transformation = from_type.make_transformation(to_type)
        result = transformation(view)
        result.path._move_or_copy(self.path_maker(**kwargs))

    def _validate_members(self, collected_paths, level):
        found_members = False
        root = pathlib.Path(self._directory_format.path)
        for path in collected_paths:
            if re.match(self.pathspec, str(path.relative_to(root))):
                if collected_paths[path]:
                    # Not a ValidationError, this just shouldn't happen.
                    raise ValueError("%r was already validated by another"
                                     " field, the pathspecs (regexes) must"
                                     " overlap." % path)
                collected_paths[path] = True
                found_members = True
                self.format(path, mode='r').validate(level)
        if not found_members:
            raise ValidationError(
                "Missing one or more files for %s: %r"
                % (self._directory_format.__class__.__name__, self.pathspec))

    @property
    def path_maker(self):
        def bound_path_maker(**kwargs):
            # Must wrap in a naive Path, otherwise an OutPath would be summoned
            # into this world, and would destroy everything in its path.
            path = (pathlib.Path(self._directory_format.path) /
                    self._path_maker(self._directory_format, **kwargs))
            # NOTE: path makers are bound to the directory format, so must be
            # provided as the first argument which will look like `self` to
            # the plugin-dev.
            path.parent.mkdir(parents=True, exist_ok=True)
            return path
        return bound_path_maker


class BoundFileCollection(BoundFile):
    def __init__(self, name, pathspec, format, directory_format, path_maker):
        super().__init__(name, pathspec, format, directory_format)
        self._path_maker = path_maker

    def view(self, view_type):
        raise NotImplementedError("Use `iter_views` instead.")

    def iter_views(self, view_type):
        # Don't want an OutPath, just a Path
        root = pathlib.Path(self._directory_format.path)
        paths = [fp for fp in sorted(root.glob('**/*'))
                 if re.match(self.pathspec, str(fp.relative_to(root)))]
        from_type = transform.ModelType.from_view_type(self.format)
        to_type = transform.ModelType.from_view_type(view_type)

        transformation = from_type.make_transformation(to_type)
        for fp in paths:
            # TODO: include capture?
            yield fp.relative_to(root), transformation(fp)


class _DirectoryMeta(type):
    def __init__(self, name, bases, dct):
        super().__init__(name, bases, dct)
        if hasattr(self, '_fields'):
            fields = self._fields.copy()
        else:
            fields = []
        for key, value in dct.items():
            if isinstance(value, File):
                # TODO: validate that the paths described by `value` are unique
                # within a DirectoryFormat
                value.name = key
                fields.append(key)

        self._fields = fields


class DirectoryFormat(FormatBase, metaclass=_DirectoryMeta):
    def validate(self, level='max'):
        _check_validation_level(level)

        if not self.path.is_dir():
            raise ValidationError("%s is not a directory." % self.path)
        collected_paths = {p: None for p in self.path.glob('**/*')
                           if not p.name.startswith('.') and
                           p.is_file()}
        for field in self._fields:
            getattr(self, field)._validate_members(collected_paths, level)

        for path, value in collected_paths.items():
            if value:
                continue
            if value is None:
                raise ValidationError("Unrecognized file (%s) for %s."
                                      % (path, self.__class__.__name__))
        if hasattr(self, '_validate_'):
            try:
                self._validate_(level)
            except ValidationError as e:
                raise ValidationError(
                    "%s is not a(n) %s:\n\n%s"
                    % (self.path, self.__class__.__name__, str(e))
                    ) from e


class SingleFileDirectoryFormatBase(DirectoryFormat):
    pass


def SingleFileDirectoryFormat(name, pathspec, format):
    # TODO: do the same hack namedtuple does so we don't mangle globals
    # (arguably the code is going to be broken if defined dynamically anyways,
    # but better to find that out later than writing in the module namespace
    # even if it isn't called module-level [which is must be!])
    df = globals()[name] = type(name, (SingleFileDirectoryFormatBase,),
                                {'file': File(pathspec, format=format)})
    return df
