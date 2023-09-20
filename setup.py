# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

setup(
    name='qiime2',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    url='https://qiime2.org',
    packages=find_packages(),
    entry_points={
        'qiime2.plugins': [
            'dummy-plugin=qiime2.core.testing.plugin:dummy_plugin',
            'other-plugin=qiime2.core.testing.plugin:other_plugin'
        ],
        'qiime2.usage_drivers': [
            'python3=qiime2.core.archive.provenance_lib:ReplayPythonUsage'
        ]
    },
    package_data={
        'qiime2.core.archive.provenance_lib.tests': ['data/**/*'],
        'qiime2.plugin.model.tests': ['data/*/*'],
        'qiime2.metadata.tests': ['data/*/*'],
        'qiime2.core.testing': ['citations.bib'],
        'qiime2.sdk.tests': ['data/*'],
        'qiime2': ['citations.bib']
    },
    zip_safe=False,
)
