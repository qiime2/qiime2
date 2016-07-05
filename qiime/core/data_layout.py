# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


class DataLayout:
    def __init__(self, name, version, validator):
        self.name = name
        self.version = version
        self._validator = validator
        self._reader_views = {}
        self._writer_views = {}

    @property
    def readers(self):
        return self._reader_views

    @property
    def writers(self):
        return self._writer_views

    def reader(self, view_type):
        def decorator(reader_function):
            self._reader_views[view_type] = reader_function
            return reader_function
        return decorator

    def writer(self, view_type):
        def decorator(writer_function):
            self._writer_views[view_type] = writer_function
            return writer_function
        return decorator

    def validate(self, data_dir):
        return self._validator(data_dir)
