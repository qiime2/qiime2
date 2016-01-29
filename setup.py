# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import glob
from setuptools import find_packages, setup

setup(
    name='qiime',
    version='2.0.1.dev0',
    packages=find_packages(),
    scripts=glob.glob('scripts/*'),
    install_requires=['click', 'python-frontmatter', 'pyyaml',
                      # Temporary dependencies:
                      'tornado'],
    package_data={'qiime': ['interface/q2d3/static/*']}
)
