# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import abc

from .base import FormatBase


class _FileFormat(FormatBase, metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def sniff(self):
        pass

    def validate(self):
        if not self.path.is_file():
            raise ValueError("%r is not a file." % self.path)
        try:
            is_member = self.sniff()
        except Exception:
            raise ValueError("Failed to sniff %r as %s. There may be a"
                             " problem with the sniffer."
                             % (self.path, self.__class__.__name__))

        if not is_member:
            raise ValueError("%r is not formatted as a %s file."
                             % (self.path, self.__class__.__name__))


class TextFileFormat(_FileFormat):
    def open(self):
        mode = 'r' if self._mode == 'r' else 'r+'
        return self.path.open(mode=mode, encoding='utf8')


class BinaryFileFormat(_FileFormat):
    def open(self):
        mode = 'rb' if self._mode == 'r' else 'r+b'
        return self.path.open(mode=mode)
