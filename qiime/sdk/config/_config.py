# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
import click
import configparser
import os

CONFIG_DIR = click.get_app_dir('QIIME 2', roaming=False, force_posix=False)
VERSION_DIR = os.path.join(CONFIG_DIR, '2.0.x')
CONFIG_FILE = os.path.join(VERSION_DIR, 'config.ini')

CONFIG = configparser.ConfigParser()
CONFIG.read(CONFIG_FILE)
