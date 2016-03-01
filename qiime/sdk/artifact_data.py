# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import abc
import tarfile
import tempfile


class ArtifactDataBase(metaclass=abc.ABCMeta):
    _root = 'data'

    @abc.abstractmethod
    def get_paths(self):
        """Return all paths to artifact data."""
        pass


class ArtifactDataReader(ArtifactDataBase):
    def __init__(self, tarfilepath):
        self._tarfile = tarfile.open(tarfilepath, mode='r')

    def get_paths(self):
        all_paths = self._tarfile.getnames()
        paths = []
        for path in all_paths:
            if os.path.dirname(path) == self._root:
                paths.append(os.path.basename(path))
        return paths

    def get_file(self, path):
        """Returns ``io.BufferedReader``"""
        filehandle = self._tarfile.extractfile(os.path.join(self._root, path))
        if filehandle is None:
            raise FileNotFoundError(
                "Filepath %r does not exist or is not a file" % path)
        return filehandle


class ArtifactDataWriter(ArtifactDataBase):
    def __init__(self):
        # TODO pass temp dir specified by user in config
        self._tempdir = tempfile.TemporaryDirectory(
            prefix='qiime2-temp-artifact-data-')
        self._tracked_files = {}

    def get_paths(self):
        # TODO this isn't ordered; does it matter?
        return list(self._tracked_files)

    def create_file(self, path):
        """Returns ``io.BufferedWriter``"""
        if path in self._tracked_files:
            raise FileExistsError("Filepath %r already exists" % path)

        filehandle = open(os.path.join(self._tempdir.name, path), mode='wb')
        self._tracked_files[path] = filehandle
        return filehandle

    def save(self, tarfilepath):
        """Returns ``tarfile.TarFile``"""
        # TODO support compression
        archive = tarfile.open(tarfilepath, mode='x')

        for filehandle in self._tracked_files.values():
            filehandle.close()

        archive.add(self._tempdir.name, arcname=self._root)
        return archive
