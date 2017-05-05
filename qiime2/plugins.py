# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import importlib.machinery

__all__ = ['available_plugins']
__path__ = []


def available_plugins():
    import qiime2.sdk
    pm = qiime2.sdk.PluginManager()
    return set('qiime2.plugins.' + s.replace('-', '_') for s in pm.plugins)


class QIIMEArtifactAPIImporter:
    def _plugin_lookup(self, plugin_name):
        import qiime2.sdk
        pm = qiime2.sdk.PluginManager()
        lookup = {s.replace('-', '_'): s for s in pm.plugins}
        if plugin_name not in lookup:
            return None
        return pm.plugins[lookup[plugin_name]]

    def find_spec(self, name, path=None, target=None):
        # Don't waste time doing anything if it's not a qiime2 plugin
        if not name.startswith('qiime2.plugins.'):
            return None

        if target is not None:
            # TODO: experiment with this to see if it is possible
            raise ImportError("Reloading the QIIME 2 Artifact API is not"
                              " currently supported.")

        # We couldn't care less about path, it is useless to us
        # (It is the __path__ of the parent module)

        fqn = name.split('.')
        plugin_details = fqn[2:]  # fqn[len(['qiime2', 'plugins']):]
        plugin_name = plugin_details[0]

        plugin = self._plugin_lookup(plugin_name)
        if plugin is None or len(plugin_details) > 2:
            return None

        if len(plugin_details) == 1:
            return self._make_spec(name, plugin)
        elif plugin_details[1] == 'visualizers':
            return self._make_spec(name, plugin, ('visualizers',))
        elif plugin_details[1] == 'methods':
            return self._make_spec(name, plugin, ('methods',))
        elif plugin_details[1] == 'actions':
            return self._make_spec(name, plugin, ('methods', 'visualizers'))
        return None

    def _make_spec(self, name, plugin, action_types=None):
        # See PEP 451 for explanation of what is happening:
        # https://www.python.org/dev/peps/pep-0451/#modulespec
        return importlib.machinery.ModuleSpec(
            name,
            loader=self,
            origin='generated QIIME 2 API',
            loader_state={'plugin': plugin, 'action_types': action_types},
            is_package=action_types is None
        )

    def create_module(self, spec):
        # Required by Python 3.6, we just need the default behavior
        return None

    def exec_module(self, module):
        spec = module.__spec__
        plugin = spec.loader_state['plugin']
        action_types = spec.loader_state['action_types']
        module.__plugin__ = plugin

        if action_types is None:
            module.methods = importlib.import_module('.methods',
                                                     package=spec.name)
            module.visualizers = importlib.import_module('.visualizers',
                                                         package=spec.name)
            module.actions = importlib.import_module('.actions',
                                                     package=spec.name)
        else:
            for action_type in action_types:
                actions = getattr(plugin, action_type)
                for key, value in actions.items():
                    setattr(module, key, value)


sys.meta_path += [QIIMEArtifactAPIImporter()]
