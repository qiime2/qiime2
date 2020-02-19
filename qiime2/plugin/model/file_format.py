# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import abc

from qiime2.core import transform
from .base import FormatBase, ValidationError, _check_validation_level


class _FileFormat(FormatBase, metaclass=abc.ABCMeta):

    def validate(self, level='max'):
        _check_validation_level(level)

        if not self.path.is_file():
            raise ValidationError("%s is not a file." % self.path)

        if hasattr(self, '_validate_'):
            try:
                self._validate_(level)
            except ValidationError as e:
                raise ValidationError(
                    "%s is not a(n) %s file:\n\n%s"
                    % (self.path, self.__class__.__name__, str(e))
                    ) from e
        # TODO: remove this branch
        elif hasattr(self, 'sniff'):
            if not self.sniff():
                raise ValidationError("%s is not a(n) %s file"
                                      % (self.path, self.__class__.__name__))

        # TODO: define an abc.abstractmethod for `validate` when sniff is
        # removed instead of this
        else:
            raise NotImplementedError("%r does not implement validate."
                                      % type(self))

    def view(self, view_type):
        from_type = transform.ModelType.from_view_type(self.__class__)
        to_type = transform.ModelType.from_view_type(view_type)

        transformation = from_type.make_transformation(to_type)
        return transformation(self)


class TextFileFormat(_FileFormat):
    def open(self):
        mode = 'r' if self._mode == 'r' else 'r+'
        # ignore BOM only when reading, do not emit BOM on write
        encoding = 'utf-8-sig' if mode == 'r' else 'utf-8'
        return self.path.open(mode=mode, encoding=encoding)


class BinaryFileFormat(_FileFormat):
    def open(self):
        mode = 'rb' if self._mode == 'r' else 'r+b'
        return self.path.open(mode=mode)
