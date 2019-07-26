# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import collections
import distutils.dir_util
import pathlib

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
        archiver = archive.Archiver.load(filepath)

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

    def save(self, filepath):
        if not filepath.endswith(self.extension):
            filepath += self.extension
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
        to_type = transform.ModelType.from_view_type(view_type)

        transformation = from_type.make_transformation(to_type,
                                                       recorder=recorder)
        result = transformation(self._archiver.data_dir)

        if view_type is qiime2.Metadata:
            result._add_source_artifact(self)

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
