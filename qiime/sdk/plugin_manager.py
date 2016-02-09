# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources


class PluginManager(object):

    def __init__(self):
        self.plugins = {}
        for plugin in pkg_resources.iter_entry_points(group='qiime.plugin'):
            self.plugins[plugin.name] = plugin.load()
