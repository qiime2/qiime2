# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import distutils.dir_util
import uuid
import glob
import os.path

import qiime.core.archiver
import qiime.core.result_base
import qiime.core.type


class Visualization(qiime.core.result_base.ResultBase):
    @classmethod
    def _assert_valid_type(cls, type_):
        if type_ is not qiime.core.type.Visualization:
            error_suffix = ""
            if qiime.core.type.is_semantic_type(type_) and type_.is_concrete():
                error_suffix = " Use qiime.sdk.Artifact with this type."
            raise TypeError(
                "A visualization requires type %r, not type %r.%s" %
                (qiime.core.type.Visualization, type_, error_suffix))

    @classmethod
    def _from_data_dir(cls, data_dir, provenance):
        viz = cls.__new__(cls)
        viz._archiver = qiime.core.archiver.Archiver(
            uuid.uuid4(), qiime.core.type.Visualization, provenance)
        # shutil.copytree doesn't allow the destination directory to exist.
        viz._archiver.save_data(data_dir, distutils.dir_util.copy_tree)
        return viz

    def get_indices(self):
        indices = glob.glob(os.path.join(self._archiver._data_dir, 'index.*'))
        return {os.path.join('data', os.path.split(fp)[1]) for fp in indices}
