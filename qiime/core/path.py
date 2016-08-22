# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import pathlib
import shutil
import tempfile
import typing


T = typing.TypeVar('T')
_ConcretePath = type(pathlib.Path())


class TempPath(_ConcretePath, typing.Generic[T]):
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
