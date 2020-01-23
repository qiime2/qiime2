# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.core.path as qpath


class FormatBase:
    def __init__(self, path=None, mode='w'):
        import qiime2.plugin.model as model
        if path is None:
            if mode != 'w':
                raise ValueError("A path must be provided when reading.")
        else:
            if mode != 'r':
                raise ValueError("A path must be omitted when writing.")

        if mode == 'w':
            self.path = qpath.OutPath(
                # TODO: parents shouldn't know about their children
                dir=isinstance(self, model.DirectoryFormat),
                prefix='q2-%s-' % self.__class__.__name__)
        else:
            self.path = qpath.InPath(path)

        self._mode = mode

    def __str__(self):
        return str(self.path)
