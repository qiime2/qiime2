# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

setup(
    name='qiime2',
    version='2017.3.0.dev',
    license='BSD-3-Clause',
    url='https://qiime2.org',
    packages=find_packages(),
    install_requires=['pyyaml', 'decorator', 'pandas', 'tzlocal',
                      'python-dateutil'],
    entry_points={
        'qiime2.plugins': [
            'dummy-plugin=qiime2.core.testing.plugin:dummy_plugin'
        ]
    },
    package_data={
        'qiime2.tests': ['data/*']
    }
)
