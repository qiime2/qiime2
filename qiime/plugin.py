# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import pkg_resources

import qiime.sdk
from qiime.core.type.type import BaseType
from qiime.core.type.primitive import Int, Str, Float

__all__ = ['Plugin', 'Type', 'Int', 'Str', 'Float']


class Plugin:
    def __init__(self, name, version, website, package):
        self.package = package
        self.name = name
        self.version = version
        self.website = website
        self.workflows = {}
        self.archive_formats = {}

    def register_workflow(self, workflow):
        fn = pkg_resources.resource_filename(self.package, workflow)
        w = qiime.sdk.Workflow.from_markdown(fn)
        self.workflows[w.id] = w

    def register_function(self, name, function, inputs, outputs, doc=""):
        # TODO where is the best place to convert outputs as a list of tuples
        # into an OrderedDict?
        outputs = collections.OrderedDict(outputs)
        w = qiime.sdk.Workflow.from_function(function, inputs, outputs, name,
                                             doc)
        self.workflows[w.id] = w

    def register_archive_format(self, name, version, validator):
        self.archive_formats[(name, version)] = \
            ArchiveFormat(name, version, validator)

    def __eq__(self, other):
        return (
            self.package == other.package and
            self.name == other.name and
            self.version == other.version and
            self.website == other.website and
            self.workflows == other.workflows and
            self.archive_formats == other.archive_formats
        )


class Type(BaseType, fields=('Artifact', 'Metadata')):
    class Artifact:
        def save(self, data, data_writer):
            pass

        def load(self, data_reader):
            pass

    class Metadata:
        def get_columns(self, data):
            pass

        def get_series(self, data, column):
            pass


class ArchiveFormat:

    def __init__(self, name, version, validator):
        # TODO: should this constructor be private since calling it
        # directly would return an unregistered format, and what
        # would you do with that? we're considering this private for
        # now, probably revisit __eq__ if that changes.
        self.name = name
        self.version = version
        self._validator = validator
        self._reader_views = {}
        self._writer_views = {}

    def __eq__(self, other):
        return (self.name == other.name and
                self.version == other.version)

    def __ne__(self, other):
        return not self == other

    def validate(self, data_reader):
        return self._validator(data_reader)

    def get_reader_views(self):
        # returns view types that this format can be read into
        return set(self._reader_views)

    def get_writer_views(self):
        # returns view types that this format can be written to
        return set(self._writer_views)

    def reader(self, view_type):
        def decorator(reader_function):
            self._reader_views[view_type] = reader_function
            return reader_function
        return decorator

    @property
    def readers(self):
        return self._reader_views

    def writer(self, view_type):
        def decorator(writer_function):
            self._writer_views[view_type] = writer_function
            return writer_function
        return decorator

    @property
    def writers(self):
        return self._writer_views
