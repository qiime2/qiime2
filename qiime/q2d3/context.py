# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


import tempfile
import glob

from qiime.sdk.system_context import SystemContext
from qiime.sdk.artifact import Artifact


class Q2D3Context(SystemContext):

    _file_extension = '.qtf'

    def __init__(self, data_dir):
        self._data_dir = data_dir
        # uuid to filepath
        data_files = glob.glob(os.path.join(
            [self._data_dir, '*%s' % _file_extension]))
        self.data = {Artifact(fp).uuid: fp for fp in data_files}

    def _output_to_filepath(self, name, output_type, workflow_template,
                            artifacts, parameters):
        with tempfile.NamedTemporaryFile(dir=self._data_dir, delete=True,
                                         suffix=_file_extension) as fh:
            return fh.name

    def _uuid_to_filepath(self, uuid):
        return self.data[uuid]
