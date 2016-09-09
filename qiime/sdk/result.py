# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import distutils.dir_util
import functools
import os
import os.path
import uuid
import pathlib

import qiime.sdk
import qiime.core.archiver as archiver
import qiime.core.type
import qiime.core.util as util
import qiime.core.transform as transform
import qiime.plugin.model as model

# Note: Result, Artifact, and Visualization classes are in this file to avoid
# circular dependencies between Result and its subclasses. Result is tightly
# coupled to Artifact and Visualization because it is a base class and a
# factory, so having the classes in the same file helps make this coupling
# explicit.


ResultMetadata = collections.namedtuple('ResultMetadata',
                                        ['uuid', 'type', 'format',
                                         'provenance'])


class Result:
    """Base class for QIIME 2 result classes (Artifact and Visualization).

    This class is not intended to be instantiated. Instead, it acts as a public
    factory and namespace for interacting with Artifacts and Visualizations in
    a generic way. It also acts as a base class for code reuse and provides an
    API shared by Artifact and Visualization.

    """

    # Subclasses must override to provide a file extension.
    extension = None

    @classmethod
    def _is_valid_type(cls, type_):
        """Subclasses should override this method."""
        return True

    @classmethod
    def peek(cls, filepath):
        uuid, type, format, provenance = archiver.Archiver.peek(filepath)
        if not cls._is_valid_type(type):
            raise TypeError(
                "Cannot peek at filepath %r because %s does not support type "
                "%r." % (filepath, cls.__name__, type))
        return ResultMetadata(uuid=uuid, type=type, format=format,
                              provenance=provenance)

    @classmethod
    def load(cls, filepath):
        """Factory for loading Artifacts and Visualizations."""
        archiver_ = archiver.Archiver.load(filepath)

        if Artifact._is_valid_type(archiver_.type):
            result = Artifact.__new__(Artifact)
        elif Visualization._is_valid_type(archiver_.type):
            result = Visualization.__new__(Visualization)
        else:
            raise TypeError(
                "Cannot load filepath %r into an Artifact or Visualization "
                "because type %r is not supported."
                % (filepath, archiver_.type))

        if type(result) is not cls and cls is not Result:
            raise TypeError(
                "Attempting to load %s with `%s.load`. Use `%s.load` instead."
                % (type(result).__name__, cls.__name__,
                   type(result).__name__))

        result._archiver = archiver_
        return result

    @classmethod
    def extract(cls, filepath, output_dir):
        # `peek` first to ensure the type is valid for the `cls` this was
        # called on.
        cls.peek(filepath)
        return archiver.Archiver.extract(filepath, output_dir)

    @property
    def type(self):
        return self._archiver.type

    @property
    def provenance(self):
        return self._archiver.provenance

    @property
    def uuid(self):
        return self._archiver.uuid

    @property
    def format(self):
        # Not memoized because it matters, just to keep the logic contained
        # in lieu of a better place to put it.
        # TODO: parse in the archiver if possible, but it will require a
        # Visualization format (currently a string), otherwise special-casing
        # is needed everywhere.
        if not hasattr(self, '__format'):
            self.__format = util.parse_format(self._archiver.format)
        return self.__format

    def __init__(self):
        raise NotImplementedError(
            "%(classname)s constructor is private, use `%(classname)s.load`, "
            "`%(classname)s.peek`, or `%(classname)s.extract`."
            % {'classname': self.__class__.__name__})

    def __new__(cls):
        result = object.__new__(cls)
        result._archiver = None
        return result

    def __repr__(self):
        return ("<%s: %r uuid: %s>"
                % (self.__class__.__name__.lower(), self.type, self.uuid))

    def __eq__(self, other):
        # Checking the UUID is mostly sufficient but requiring an exact type
        # match makes it safer in case `other` is a subclass or a completely
        # different type that happens to have a `.uuid` property. We want to
        # ensure (as best as we can) that the UUIDs we are comparing are linked
        # to the same type of QIIME 2 object.
        return (
            type(self) == type(other) and
            self.uuid == other.uuid
        )

    def __ne__(self, other):
        return not (self == other)

    def _orphan(self, pid):
        self._archiver.orphan(pid)

    def save(self, filepath):
        if not filepath.endswith(self.extension):
            filepath += self.extension
        self._archiver.save(filepath)
        return filepath


class Artifact(Result):
    extension = '.qza'

    @classmethod
    def _is_valid_type(cls, type_):
        if qiime.core.type.is_semantic_type(type_) and type_.is_concrete():
            return True
        else:
            return False

    @classmethod
    def import_data(cls, type, view, view_type=None):
        type_, type = type, __builtins__['type']
        if isinstance(type_, str):
            type_ = qiime.sdk.parse_type(type_)

        if isinstance(view_type, str):
            view_type = qiime.core.util.parse_format(view_type)

        if view_type is None:
            if type(view) is str or isinstance(view, pathlib.PurePath):
                pm = qiime.sdk.PluginManager()
                output_dir_fmt = pm.get_directory_format(type_)
                if pathlib.Path(view).is_file() and issubclass(
                        output_dir_fmt,
                        model.SingleFileDirectoryFormatBase):
                    view_type = output_dir_fmt.file.format
                else:
                    view_type = output_dir_fmt
            else:
                view_type = type(view)

        provenance = ("This artifact was generated by importing data from an "
                      "external source.")
        return cls._from_view(type_, view, view_type, provenance)

    @classmethod
    def _from_view(cls, type, view, view_type, provenance):
        if isinstance(type, str):
            type = qiime.sdk.parse_type(type)

        if not cls._is_valid_type(type):
            raise TypeError(
                "An artifact requires a concrete semantic type, not type %r."
                % type)

        pm = qiime.sdk.PluginManager()
        output_dir_fmt = pm.get_directory_format(type)

        if view_type is None:
            # lookup default format for the type
            view_type = output_dir_fmt

        from_type = transform.ModelType.from_view_type(view_type)
        to_type = transform.ModelType.from_view_type(output_dir_fmt)

        transformation = from_type.make_transformation(to_type)
        result = transformation(view)

        artifact = cls.__new__(cls)
        artifact._archiver = archiver.Archiver(
            uuid.uuid4(), type, output_dir_fmt.__name__, provenance,
            data_initializer=result.path._move_or_copy)
        return artifact

    def view(self, view_type):
        from_type = transform.ModelType.from_view_type(self.format)
        to_type = transform.ModelType.from_view_type(view_type)
        transformation = from_type.make_transformation(to_type)
        result = self._archiver.load_data(transformation)

        to_type.set_user_owned(result, True)
        return result


class Visualization(Result):
    extension = '.qzv'

    @classmethod
    def _is_valid_type(cls, type_):
        if type_ == qiime.core.type.Visualization:
            return True
        else:
            return False

    @classmethod
    def _from_data_dir(cls, data_dir, provenance):
        # shutil.copytree doesn't allow the destination directory to exist.
        data_initializer = functools.partial(distutils.dir_util.copy_tree,
                                             data_dir)
        viz = cls.__new__(cls)
        viz._archiver = archiver.Archiver(
            uuid.uuid4(), qiime.core.type.Visualization, 'Visualization',
            provenance, data_initializer=data_initializer)
        return viz

    def get_index_paths(self, relative=True):
        result = {}
        for relpath, abspath in self._archiver.get_data_paths(recursive=False):
            if relpath.startswith('index.'):
                relpath = os.path.join(self._archiver.DATA_DIRNAME, relpath)
                ext = os.path.splitext(relpath)[1][1:]
                if ext in result:
                    raise ValueError(
                        "Multiple index files identified with %s "
                        "extension (%s, %s). This is currently "
                        "unsupported." %
                        (ext, result[ext], relpath))
                else:
                    result[ext] = relpath if relative else abspath
        return result
