# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import os

from qiime.sdk import Artifact


class Q2D3Context:

    _file_extension = '.qtf'

    def __init__(self, data_dir, output_names=None):
        self._data_dir = data_dir
        # uuid to filepath
        data_files = glob.glob(os.path.join(
            self._data_dir, '*%s' % self._file_extension))
        self.data = {Artifact(fp).uuid: fp for fp in data_files}
        self.names = {
            Artifact(fp).uuid: os.path.splitext(os.path.split(fp)[1])[0]
            for fp in data_files
        }

        self.output_names = output_names
