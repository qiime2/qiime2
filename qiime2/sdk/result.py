# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import warnings
import tempfile
import collections
import distutils.dir_util
import pathlib
from typing import Union, get_args, get_origin

import qiime2.metadata
import qiime2.plugin
import qiime2.sdk
import qiime2.core.type
import qiime2.core.transform as transform
import qiime2.core.archive as archive
import qiime2.plugin.model as model
import qiime2.core.util as util
import qiime2.core.exceptions as exceptions

# Note: Result, Artifact, and Visualization classes are in this file to avoid
# circular dependencies between Result and its subclasses. Result is tightly
# coupled to Artifact and Visualization because it is a base class and a
# factory, so having the classes in the same file helps make this coupling
# explicit.


ResultMetadata = collections.namedtuple('ResultMetadata',
                                        ['uuid', 'type', 'format'])


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
        return ResultMetadata(*archive.Archiver.peek(filepath))

    @classmethod
    def extract(cls, filepath, output_dir):
        """Unzip contents of Artifacts and Visualizations."""
        return archive.Archiver.extract(filepath, output_dir)

    @classmethod
    def load(cls, filepath):
        """Factory for loading Artifacts and Visualizations."""
        from qiime2.core.cache import get_cache

        # Check if the data is already in the cache (if the uuid is in
        # cache.data) and load it from the cache if it is. Avoids unzipping the
        # qza again if we already have it.
        cache = get_cache()
        peek = cls.peek(filepath)
        archiver = cache._load_uuid(peek.uuid)

        if not archiver:
            try:
                archiver = archive.Archiver.load(filepath)
            except OSError as e:
                if e.errno == 28:
                    temp = tempfile.tempdir
                    raise ValueError(f'There was not enough space left on '
                                     f'{temp!r} to extract the artifact '
                                     f'{filepath!r}. (Try setting $TMPDIR to '
                                     'a directory with more space, or '
                                     f'increasing the size of {temp!r})')
                else:
                    raise e

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
    def _from_archiver(cls, archiver):
        if Artifact._is_valid_type(archiver.type):
            result = Artifact.__new__(Artifact)
        elif Visualization._is_valid_type(archiver.type):
            result = Visualization.__new__(Visualization)
        else:
            raise TypeError(
                "Cannot load filepath %r into an Artifact or Visualization "
                "because type %r is not supported."
                % (archiver.path, archiver.type))

        if type(result) is not cls and cls is not Result:
            raise TypeError(
                "Attempting to load %s with `%s.load`. Use `%s.load` instead."
                % (type(result).__name__, cls.__name__,
                   type(result).__name__))

        result._archiver = archiver
        return result

    @property
    def type(self):
        return self._archiver.type

    @property
    def uuid(self):
        return self._archiver.uuid

    @property
    def format(self):
        return self._archiver.format

    @property
    def citations(self):
        return self._archiver.citations

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

    def __hash__(self):
        return hash(self.uuid)

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

    def export_data(self, output_dir):
        distutils.dir_util.copy_tree(
            str(self._archiver.data_dir), str(output_dir))
        # Return None for now, although future implementations that include
        # format tranformations may return the invoked transformers
        return None

    @property
    def _destructor(self):
        return self._archiver._destructor

    def save(self, filepath, ext=None):
        """Save to a file.

        Parameters
        ----------
        filepath : str
            Path to save file at.

        extension : str
            Preferred file extension (.qza, .qzv, .txt, etc).
            If no preferred extension input is included,
            Artifact extension will default to .qza and
            Visualization extension will default to .qzv.
            Including a period in the extension is
            optional, and any additional periods delimiting
            the filepath and the extension will be reduced
            to a single period.

        Returns
        -------
        str
            Filepath and extension (if provided) that the
            file was saved to.

        See Also
        --------
        load

        """
        if ext is None:
            ext = self.extension

        # This accounts for edge cases in the filename extension
        # and ensures that there is only a single period in the ext.
        filepath = filepath.rstrip('.')
        ext = '.' + ext.lstrip('.')

        if not filepath.endswith(ext):
            filepath += ext

        self._archiver.save(filepath)
        return filepath

    def _alias(self, provenance_capture):
        def clone_original(into):
            # directory is empty, this function is meant to fix that, so we
            # can rmdir so that copytree is happy
            into.rmdir()
            shutil.copytree(str(self._archiver.data_dir), str(into),
                            copy_function=os.link)  # Use hardlinks

        cls = type(self)
        alias = cls.__new__(cls)
        alias._archiver = archive.Archiver.from_data(
            self.type, self.format, clone_original, provenance_capture)
        return alias

    def validate(self, level=NotImplemented):
        diff = self._archiver.validate_checksums()
        if diff.changed or diff.added or diff.removed:
            error = ""

            if diff.added:
                error += "Unrecognized files:\n"
            for key in diff.added:
                error += "  - %r\n" % key
            if diff.removed:
                error += "Missing files:\n"
            for key in diff.removed:
                error += "  - %r\n" % key
            if diff.changed:
                error += "Changed files:\n"
            for (key, (exp, obs)) in diff.changed.items():
                error += "  - %r: %s -> %s\n" % (key, exp, obs)

            raise exceptions.ValidationError(error)

    def result(self):
        """ Noop to provide standardized interface with ProxyResult.
        """
        return self


class Artifact(Result):
    extension = '.qza'

    @classmethod
    def _is_valid_type(cls, type_):
        if qiime2.core.type.is_semantic_type(type_) and type_.is_concrete():
            return True
        else:
            return False

    @classmethod
    def import_data(cls, type, view, view_type=None):
        type_, type = type, __builtins__['type']

        is_format = False
        if isinstance(type_, str):
            type_ = qiime2.sdk.parse_type(type_)

        if isinstance(view_type, str):
            view_type = qiime2.sdk.parse_format(view_type)
            is_format = True

        if view_type is None:
            if type(view) is str or isinstance(view, pathlib.PurePath):
                is_format = True
                pm = qiime2.sdk.PluginManager()
                output_dir_fmt = pm.get_directory_format(type_)
                if pathlib.Path(view).is_file():
                    if not issubclass(output_dir_fmt,
                                      model.SingleFileDirectoryFormatBase):
                        raise qiime2.plugin.ValidationError(
                            "Importing %r requires a directory, not %s"
                            % (output_dir_fmt.__name__, view))
                    view_type = output_dir_fmt.file.format
                else:
                    view_type = output_dir_fmt
            else:
                view_type = type(view)

        format_ = None
        md5sums = None
        if is_format:
            path = pathlib.Path(view)
            if path.is_file():
                md5sums = {path.name: util.md5sum(path)}
            elif path.is_dir():
                md5sums = util.md5sum_directory(path)
            else:
                raise qiime2.plugin.ValidationError(
                    "Path '%s' does not exist." % path)
            format_ = view_type

        provenance_capture = archive.ImportProvenanceCapture(format_, md5sums)
        return cls._from_view(type_, view, view_type, provenance_capture,
                              validate_level='max')

    @classmethod
    def _from_view(cls, type, view, view_type, provenance_capture,
                   validate_level='min'):
        type_raw = type
        if isinstance(type, str):
            type = qiime2.sdk.parse_type(type)

        if not cls._is_valid_type(type):
            raise TypeError(
                "An artifact requires a concrete semantic type, not type %r."
                % type)

        pm = qiime2.sdk.PluginManager()
        output_dir_fmt = pm.get_directory_format(type)

        if view_type is None:
            # lookup default format for the type
            view_type = output_dir_fmt

        from_type = transform.ModelType.from_view_type(view_type)
        to_type = transform.ModelType.from_view_type(output_dir_fmt)

        recorder = provenance_capture.transformation_recorder('return')
        transformation = from_type.make_transformation(to_type,
                                                       recorder=recorder)
        result = transformation(view, validate_level)

        if type_raw in pm.validators:
            validation_object = pm.validators[type]
            validation_object(data=result, level=validate_level)

        artifact = cls.__new__(cls)
        artifact._archiver = archive.Archiver.from_data(
            type, output_dir_fmt,
            data_initializer=result.path._move_or_copy,
            provenance_capture=provenance_capture)
        return artifact

    def view(self, view_type):
        return self._view(view_type)

    def _view(self, view_type, recorder=None):
        if view_type is qiime2.Metadata and not self.has_metadata():
            raise TypeError(
                "Artifact %r cannot be viewed as QIIME 2 Metadata." % self)

        from_type = transform.ModelType.from_view_type(self.format)

        if isinstance(get_origin(view_type), type(Union)):
            transformation = None
            for arg in get_args(view_type):
                to_type = transform.ModelType.from_view_type(arg)
                try:
                    transformation = from_type.make_transformation(
                        to_type, recorder=recorder)
                    if transformation:
                        break
                except Exception as e:
                    if str(e).startswith("No transformation from"):
                        continue
                    else:
                        raise e
            if not transformation:
                raise Exception(
                    "No transformation into either of %s was found" %
                    ", ".join([str(x) for x in view_type.__args__])
                )
        else:
            to_type = transform.ModelType.from_view_type(view_type)
            transformation = from_type.make_transformation(to_type,
                                                           recorder=recorder)
        result = transformation(self._archiver.data_dir)

        if view_type is qiime2.Metadata:
            result._add_artifacts([self])

        to_type.set_user_owned(result, True)
        return result

    def has_metadata(self):
        """ Checks for metadata within an artifact

        Returns
        -------
        bool
           True if the artifact has metadata (i.e. can be viewed as
           ``qiime2.Metadata``), False otherwise.

        """
        from_type = transform.ModelType.from_view_type(self.format)
        to_type = transform.ModelType.from_view_type(qiime2.Metadata)
        return from_type.has_transformation(to_type)

    def validate(self, level='max'):
        """ Validates the data contents of an artifact

        Raises
        ------
        ValidationError
            If the artifact is invalid at the specified level of validation.
        """
        super().validate()

        self.format.validate(self.view(self.format), level)


class Visualization(Result):
    extension = '.qzv'

    @classmethod
    def _is_valid_type(cls, type_):
        return type_ == qiime2.core.type.Visualization

    @classmethod
    def _from_data_dir(cls, data_dir, provenance_capture):
        # shutil.copytree doesn't allow the destination directory to exist.
        def data_initializer(destination):
            return distutils.dir_util.copy_tree(
                str(data_dir), str(destination))

        viz = cls.__new__(cls)
        viz._archiver = archive.Archiver.from_data(
            qiime2.core.type.Visualization, None,
            data_initializer=data_initializer,
            provenance_capture=provenance_capture)
        return viz

    def get_index_paths(self, relative=True):
        result = {}
        for abspath in self._archiver.data_dir.iterdir():
            data_path = str(abspath.relative_to(self._archiver.data_dir))
            if data_path.startswith('index.'):
                relpath = abspath.relative_to(self._archiver.root_dir)
                ext = relpath.suffix[1:]
                if ext in result:
                    raise ValueError(
                        "Multiple index files identified with %s "
                        "extension (%s, %s). This is currently "
                        "unsupported." %
                        (ext, result[ext], relpath))
                else:
                    result[ext] = str(relpath) if relative else str(abspath)
        return result

    def _repr_html_(self):
        from qiime2.jupyter import make_html
        return make_html(str(self._archiver.path))


class ResultCollection:
    @classmethod
    def load(cls, directory):
        """ Determines how to load a Collection of QIIME 2 Artifacts in a
            directory and dispatches to helpers
        """
        if not os.path.isdir(directory):
            raise ValueError(
                f"Given filepath '{directory}' is not a directory")

        order_fp = os.path.join(directory, '.order')

        if os.path.isfile(order_fp):
            collection = cls._load_ordered(directory, order_fp)
        else:
            warnings.warn(f"The directory '{directory}' does not contain a "
                          ".order file. The files will be read into the "
                          "collection in the order the filesystem provides "
                          "them in.")
            collection = cls._load_unordered(directory)

        return collection

    @classmethod
    def _load_ordered(cls, directory, order_fp):
        collection = cls()

        with open(order_fp, 'r') as order_fh:
            for result_name in order_fh.read().splitlines():
                result_fp = cls._get_result_fp(directory, result_name)
                collection[result_name] = Result.load(result_fp)

        return collection

    @classmethod
    def _load_unordered(cls, directory):
        collection = cls()

        for result in os.listdir(directory):
            result_fp = os.path.join(directory, result)
            result_name = result.rstrip('.qza')
            result_name = result_name.rstrip('.qzv')

            collection[result_name] = Result.load(result_fp)

        return collection

    @classmethod
    def _get_result_fp(cls, directory, result_name):
        result_fp = os.path.join(directory, result_name)

        # Check if thing in .order file exists and if not try it with .qza at
        # the end and if not try it with .qzv at the end
        if not os.path.isfile(result_fp):
            result_fp += '.qza'

            if not os.path.isfile(result_fp):
                # Get rid of the trailing .qza before adding .qzv
                result_fp = result_fp[:-4]
                result_fp += '.qzv'

                if not os.path.isfile(result_fp):
                    raise ValueError(
                        f"The Result '{result_name}' is referenced in the "
                        "order file but does not exist in the directory.")

        return result_fp

    def __init__(self, collection=None):
        if collection is None:
            self.collection = {}
        elif isinstance(collection, dict):
            self.collection = collection
        else:
            self.collection = {k: v for k, v in enumerate(collection)}

    def __contains__(self, item):
        return item in self.collection

    def __eq__(self, other):
        if isinstance(other, dict):
            return self.collection == other
        elif isinstance(other, ResultCollection):
            return self.collection == other.collection
        else:
            raise TypeError(f'Equality between {type(other)} and '
                            'ResultCollection is undefined.')

    def __len__(self):
        return len(self.collection)

    def __iter__(self):
        yield self.collection.__iter__()

    def __setitem__(self, key, item):
        self.collection[key] = item

    def __getitem__(self, key):
        return self.collection[key]

    def __repr__(self):
        return f"<{self.__class__.__name__.lower()}: {self.type}>"

    @property
    def type(self):
        inner_type = qiime2.core.type.grammar.UnionExp(
            v.type for v in self.collection.values()).normalize()

        return qiime2.core.type.Collection[inner_type]

    def save(self, directory):
        """Saves a collection of QIIME 2 Results into a given directory with
           an order file.

           NOTE: The directory given must not exist
        """
        if os.path.exists(directory):
            raise ValueError(f"The given directory '{directory}' already "
                             "exists. A new directory must be given to save "
                             "the collection to.")

        os.makedirs(directory)

        with open(os.path.join(directory, '.order'), 'w') as fh:
            for name, result in self.collection.items():
                result_fp = os.path.join(directory, name)
                result.save(result_fp)
                fh.write(f'{name}\n')

        # Do this to give us a unified API with Result.save
        return directory

    def keys(self):
        return self.collection.keys()

    def values(self):
        return self.collection.values()

    def items(self):
        return self.collection.items()

    def validate(self, view, level=None):
        for result in self.values():
            result.validate(view, level)

    def result(self):
        """ Noop to provide standardized interface with ProxyResultCollection.
        """
        return self
