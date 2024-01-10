# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
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

        from qiime2.core.cache import get_cache
        import os
        cache = get_cache()

        if mode == 'w':
            # print("HERE IN WRITE")
            self.path = qpath.OutPath(
                # TODO: parents shouldn't know about their children
                dir=isinstance(self, model.DirectoryFormat),
                prefix=os.path.join(cache.process_pool.path,
                                    'q2-%s-' % self.__class__.__name__))
        else:
            # print("HERE IN READ")
            # raise ValueError(path)
            # print(path)
            # from qiime2.core.cache import get_cache
            # import os
            # cache = get_cache()
            # path = cache.process_pool._allocate(os.path.basename(path))
            # raise ValueError(path)
            # path = os.path.join(path, 'mount')
            # print(f"PATH IN WRITE {path}")
            # time.sleep(600)
            # raise ValueError("ERROR")
            self.path = qpath.InPath(path)
            # raise ValueError(self.path)

        self._mode = mode

    def __str__(self):
        return str(self.path)
