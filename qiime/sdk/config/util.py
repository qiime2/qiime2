# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
import os
import shutil

from ._config import VERSION_DIR

__BOILERPLATE_DIR = os.path.join(os.path.dirname(__file__), 'boilerplate')


def reset_config():
    try:
        shutil.rmtree(VERSION_DIR)
    except FileNotFoundError:
        pass
    shutil.copytree(__BOILERPLATE_DIR, VERSION_DIR)


def config_location():
    print(os.path.join(VERSION_DIR, 'config.ini'))
