# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import distutils.dir_util
import functools
import uuid
import os

import qiime.core.archiver
import qiime.core.result_base
import qiime.core.type


class Visualization(qiime.core.result_base.ResultBase):
    @classmethod
    def _assert_valid_type(cls, type_):
        if type_ != qiime.core.type.Visualization:
            error_suffix = ""
            if qiime.core.type.is_semantic_type(type_) and type_.is_concrete():
                error_suffix = " Use qiime.sdk.Artifact with this type."
            raise TypeError(
                "A visualization requires type %r, not type %r.%s" %
                (qiime.core.type.Visualization, type_, error_suffix))

    @classmethod
    def _from_data_dir(cls, data_dir, provenance):
        # shutil.copytree doesn't allow the destination directory to exist.
        data_initializer = functools.partial(distutils.dir_util.copy_tree,
                                             data_dir)
        viz = cls.__new__(cls)
        viz._archiver = qiime.core.archiver.Archiver(
            uuid.uuid4(), qiime.core.type.Visualization, provenance,
            data_initializer=data_initializer)
        return viz

    def get_index_paths(self, relative=True):
        result = {}
        for relpath, abspath in self._archiver.get_data_paths(recursive=False):
            if relpath.startswith('index.'):
                relpath = os.path.join(self._archiver.DATA_DIRNAME, relpath)
                ext = os.path.splitext(relpath)[1][1:]
                if ext in result:
                    # TODO: this should [additionally] be handled
                    # elsewhere, probably on population of self._data_dir.
                    raise ValueError(
                        "Multiple index files identified with %s "
                        "extension (%s, %s). This is currently "
                        "unsupported." %
                        (ext, result[ext], relpath))
                else:
                    result[ext] = relpath if relative else abspath
        return result
