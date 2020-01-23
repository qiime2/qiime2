# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pathlib

import tornado.web as web
from notebook.base.handlers import IPythonHandler

import qiime2.core.archive.archiver as archiver


class _ArchiveCheck(archiver._Archive):
    """This is only what is needed to verify a path is an archive"""
    # TODO: make this part of the archiver API at some point
    def open(self, relpath):
        abspath = os.path.join(str(self.path), str(self.uuid), relpath)
        return open(abspath, 'r')

    def relative_iterdir(self, relpath='.'):
        for p in pathlib.Path(self.path).iterdir():
            yield str(p.relative_to(self.path))


class QIIME2RedirectHandler(IPythonHandler):
    """Add a location to location_store for later retrieval"""
    def initialize(self, result_store):
        self.result_store = result_store

    def get(self):
        location = self.get_query_argument('location')
        if not os.path.exists(location):
            # Client DOM should explain that the user should re-run the cell
            self.send_error(409)  # Conflict
            return
        # is it actually a QIIME 2 result, or a random part of the filesystem
        archive = _ArchiveCheck(pathlib.Path(location))
        self.result_store[archive.uuid] = os.path.join(
            location, str(archive.uuid), 'data')

        self.redirect('view/%s/' % archive.uuid)


class QIIME2ResultHandler(web.StaticFileHandler):
    def initialize(self, path, default_filename):
        super().initialize(path, default_filename)
        self.result_store = path  # path is actually result_store

    @classmethod
    def get_absolute_path(cls, root, path):
        uuid, path = path.split('/', 1)
        root = root[uuid]
        # This is janky, but validate_absolute_path is the only thing
        # that will use this data, so it can know to unpack the tuple again
        return (super().get_absolute_path(root, path), uuid)

    def validate_absolute_path(self, root, abspath_uuid):
        absolute_path, uuid = abspath_uuid
        root = self.result_store[uuid]
        return super().validate_absolute_path(root, absolute_path)
