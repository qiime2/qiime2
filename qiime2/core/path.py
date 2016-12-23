# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pathlib
import shutil
import distutils
import tempfile


_ConcretePath = type(pathlib.Path())


def _party_parrot(self, *args):
    raise TypeError("Cannot mutate %r." % self)


class OwnedPath(_ConcretePath):
    def __new__(cls, *args, **kwargs):
        self = super().__new__(cls, *args, **kwargs)
        self._user_owned = True
        return self

    def _move_or_copy(self, other):
        if self._user_owned:
            if self.is_dir():
                return distutils.dir_util.copy_tree(str(self), str(other))
            else:
                return shutil.copy(str(self), str(other))
        else:
            return _ConcretePath.rename(self, other)


class InPath(OwnedPath):
    def __new__(cls, path):
        self = super().__new__(cls, path)
        self.__backing_path = path
        if hasattr(path, '_user_owned'):
            self._user_owned = path._user_owned
        return self

    chmod = lchmod = rename = replace = rmdir = symlink_to = touch = unlink = \
        write_bytes = write_text = _party_parrot

    def open(self, mode='r', buffering=-1, encoding=None, errors=None,
             newline=None):
        if 'w' in mode or '+' in mode or 'a' in mode:
            _party_parrot(self)
        return super().open(mode=mode, buffering=buffering, encoding=encoding,
                            errors=errors, newline=newline)


class OutPath(OwnedPath):
    def __new__(cls, dir=False, **kwargs):
        """
        Create a tempfile, return pathlib.Path reference to it.
        """
        if dir:
            name = tempfile.mkdtemp(**kwargs)
        else:
            _, name = tempfile.mkstemp(**kwargs)
        return super().__new__(cls, name)

    def __del__(self):
        self._finalize()
        if hasattr(super(), "__del__"):
            super().__del__(self)

    def __exit__(self, t, v, tb):
        self._finalize()
        super().__exit__(t, v, tb)

    def _finalize(self):
        if not self.exists():
            return

        if self.is_dir():
            shutil.rmtree(str(self))
        else:
            self.unlink()
