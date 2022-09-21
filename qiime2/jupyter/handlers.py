# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pathlib

import tornado.web as web
from notebook.base.handlers import IPythonHandler

from qiime2.core.archive.archiver import ArchiveCheck


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
        archive = ArchiveCheck(pathlib.Path(location))
        self.result_store[archive.uuid] = os.path.join(location, 'data')

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
