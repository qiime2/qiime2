# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
import tempfile

from ._config import CONFIG


def gettempdir():
    paths = CONFIG['Paths'] if 'Paths' in CONFIG else {}
    return paths.get('temp_dir', tempfile.gettempdir())


def mkdtemp(suffix=None):
    return tempfile.mkdtemp(suffix, prefix='q2-', dir=gettempdir())


class TemporaryDirectory(tempfile.TemporaryDirectory):

    def __init__(self, suffix=None):
        super().__init__(suffix=suffix, prefix='q2-', dir=gettempdir())
