# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import qiime.sdk


def get_dummy_plugin():
    plugin_manager = qiime.sdk.PluginManager()
    if 'dummy-plugin' not in plugin_manager.plugins:
        raise RuntimeError(
            "When running QIIME 2 unit tests, the QIIMETEST environment "
            "variable must be defined so that plugins required by unit tests "
            "are loaded. The value of the QIIMETEST environment variable can "
            "be anything. Example command: QIIMETEST=1 nosetests")
    return plugin_manager.plugins['dummy-plugin']
