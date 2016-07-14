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
import shutil
import uuid

import qiime.sdk
import qiime.core.archiver
import qiime.core.type

# Note: Result, Artifact, and Visualization classes are in this file to avoid
# circular dependencies between Result and its subclasses. Result is tightly
# coupled to Artifact and Visualization because it is a base class and a
# factory, so having the classes in the same file helps make this coupling
# explicit.


ResultMetadata = collections.namedtuple('ResultMetadata',
                                        ['uuid', 'type', 'provenance'])


class Result:
    """Base class for QIIME 2 result classes (Artifact and Visualization).

    This class is not intended to be instantiated. Instead, it acts as a public
    factory and namespace for interacting with Artifacts and Visualizations in
    a generic way. It also acts as a base class for code reuse and provides an
    API shared by Artifact and Visualization.

    """
    @classmethod
    def _is_valid_type(cls, type_):
        """Subclasses should override this method."""
        return True

    @classmethod
    def peek(cls, filepath):
        uuid, type, provenance = qiime.core.archiver.Archiver.peek(filepath)
        if not cls._is_valid_type(type):
            raise TypeError(
                "Cannot peek at filepath %r because %s does not support type "
                "%r." % (filepath, cls.__name__, type))
        return ResultMetadata(uuid=uuid, type=type, provenance=provenance)

    @classmethod
    def load(cls, filepath):
        """Factory for loading Artifacts and Visualizations."""
        archiver = qiime.core.archiver.Archiver.load(filepath)

        if Artifact._is_valid_type(archiver.type):
            result = Artifact.__new__(Artifact)
        elif Visualization._is_valid_type(archiver.type):
            result = Visualization.__new__(Visualization)
        else:
            raise TypeError(
                "Cannot load filepath %r into an Artifact or Visualization "
                "because type %r is not supported."
                % (filepath, archiver.type))

        if type(result) is not cls and cls is not Result:
            raise TypeError(
                "Attempting to load %s with `%s.load`. Use `%s.load` instead."
                % (type(result).__name__, cls.__name__,
                   type(result).__name__))

        result._archiver = archiver
        return result

    @classmethod
    def extract(cls, filepath, output_dir):
        # `peek` first to ensure the type is valid for the `cls` this was
        # called on.
        cls.peek(filepath)
        return qiime.core.archiver.Archiver.extract(filepath, output_dir)

    @property
    def type(self):
        return self._archiver.type

    @property
    def provenance(self):
        return self._archiver.provenance

    @property
    def uuid(self):
        return self._archiver.uuid

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

    def _orphan(self, pid):
        self._archiver.orphan(pid)

    def save(self, filepath):
        self._archiver.save(filepath)


class Artifact(Result):
    @classmethod
    def _is_valid_type(cls, type_):
        if qiime.core.type.is_semantic_type(type_) and type_.is_concrete():
            return True
        else:
            return False

    @classmethod
    def import_data(cls, type, path):
        # Note: `import` is a better name for this method but it's a reserved
        # word in Python and can't be used as a function/method name. We use
        # `import_data` because this method will accept views (i.e. Python
        # objects) in the future in addition to files/directories. These all
        # represent an artifact's data and it is this method's responsibility
        # to transform that data into a relevant data layout.
        if isinstance(type, str):
            type = qiime.sdk.parse_type(type)
        if not cls._is_valid_type(type):
            raise TypeError(
                "An artifact requires a concrete semantic type, not type %r."
                % type)

        data_layout = cls._get_data_layout(type)
        data_layout.validate(path)

        if os.path.isfile(path):
            relpath = next(iter(data_layout.files))

            def data_initializer(data_dir):
                # TODO support data layouts with a single file in a nested
                # directory (currently this will fail if a data layout like
                # this exists).
                shutil.copyfile(path, os.path.join(data_dir, relpath))
        else:
            # shutil.copytree doesn't allow the destination directory to exist.
            data_initializer = functools.partial(
                distutils.dir_util.copy_tree, path)

        artifact = cls.__new__(cls)
        # TODO import provenance is stubbed for now as a string
        provenance = ("This artifact was generated by importing data from an "
                      "external source.")
        artifact._archiver = qiime.core.archiver.Archiver(
            uuid.uuid4(), type, provenance, data_initializer=data_initializer)
        return artifact

    @classmethod
    def _from_view(cls, view, type_, provenance):
        """

        Parameters
        ----------
        view : Python object
            View to serialize.
        type_ : qiime.plugin.Type
            Concrete semantic type of the artifact.
        provenance : qiime.sdk.Provenance
            Artifact provenance.

        """
        if not cls._is_valid_type(type_):
            raise TypeError(
                "An artifact requires a concrete semantic type, not type %r."
                % type_)
        if provenance is not None and not isinstance(provenance,
                                                     qiime.sdk.Provenance):
            raise TypeError(
                "`provenance` must be None or an instance of "
                "qiime.sdk.Provenance.")

        data_layout = cls._get_data_layout(type_)
        view_type = type(view)
        # TODO better error handling for when `view` cannot be written to
        # `type_` data layout.
        writer = data_layout.writers[view_type]
        writer = functools.partial(writer, view)

        artifact = cls.__new__(cls)
        artifact._archiver = qiime.core.archiver.Archiver(
            uuid.uuid4(), type_, provenance, data_initializer=writer)
        return artifact

    @classmethod
    def _get_data_layout(cls, type_):
        pm = qiime.sdk.PluginManager()

        data_layout = None
        for semantic_type, datalayout in \
                pm.semantic_type_to_data_layouts.items():
            if type_ <= semantic_type:
                data_layout = datalayout
                break

        if data_layout is None:
            raise TypeError(
                "Artifact semantic type %r does not have a compatible data "
                "layout." % type_)

        return data_layout

    def view(self, view_type):
        data_layout = self._get_data_layout(self.type)
        reader = data_layout.readers[view_type]
        return self._archiver.load_data(reader)


class Visualization(Result):
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
        viz._archiver = qiime.core.archiver.Archiver(
            uuid.uuid4(), qiime.core.type.Visualization, provenance,
            data_initializer=data_initializer)
        return viz

    def get_index_paths(self, relative=True):
        result = {}
        for relpath, abspath in self._archiver.get_data_paths(recursive=False):
            if relpath.startswith('index.'):
                relpath = os.path.join(self._archiver.DATA_DIRNAME, relpath)
                ext = os.path.splitext(relpath)[1][1:]
                if ext in result:
                    # TODO: this should [additionally] be handled
                    # elsewhere, probably on population of self._data_dir.
                    raise ValueError(
                        "Multiple index files identified with %s "
                        "extension (%s, %s). This is currently "
                        "unsupported." %
                        (ext, result[ext], relpath))
                else:
                    result[ext] = relpath if relative else abspath
        return result
