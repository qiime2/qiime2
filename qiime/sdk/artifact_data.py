# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tarfile
import tempfile


# TODO make sure files are closed appropriately
# TODO support subdirectories?
class ArtifactDataBase:
    _root = 'data'


class ArtifactDataReader(ArtifactDataBase):
    def __init__(self, tar):
        self._tar = tar

    def get_paths(self):
        """Return all paths to artifact data."""
        all_paths = self._tar.getnames()
        paths = []
        for path in all_paths:
            if os.path.dirname(path) == self._root:
                paths.append(os.path.basename(path))
        return paths

    def get_file(self, path):
        """Returns ``io.BufferedReader``"""
        try:
            filehandle = self._tar.extractfile(os.path.join(self._root, path))
        except KeyError:
            raise FileNotFoundError("Filepath %r does not exist" % path)

        if filehandle is None:
            raise FileNotFoundError(
                "Filepath %r is not a file" % path)
        return filehandle


class ArtifactDataWriter(ArtifactDataBase):
    def __init__(self):
        # TODO pass temp dir specified by user in config
        self._tempdir = tempfile.TemporaryDirectory(
            prefix='qiime2-temp-artifact-data-')
        self._tracked_files = {}

    def create_file(self, path):
        """Returns ``io.BufferedWriter``"""
        if path in self._tracked_files:
            raise FileExistsError("Filepath %r already exists" % path)

        filehandle = open(os.path.join(self._tempdir.name, path), mode='wb')
        self._tracked_files[path] = filehandle
        return filehandle

    def save(self, tar):
        for filehandle in self._tracked_files.values():
            filehandle.close()

        # TODO use `filter` parameter to only add files we know about
        # TODO set file metadata appropriately (e.g., owner, permissions)
        tar.add(self._tempdir.name, arcname=self._root)
