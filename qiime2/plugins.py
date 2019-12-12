# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import importlib.machinery

from qiime2.sdk import usage

__all__ = ['available_plugins', 'ArtifactAPIUsage']
__path__ = []


def available_plugins():
    import qiime2.sdk
    pm = qiime2.sdk.PluginManager()
    return set('qiime2.plugins.' + s.replace('-', '_') for s in pm.plugins)


class ArtifactAPIUsage(usage.Usage):
    def __init__(self):
        super().__init__()
        self._imports = set()
        self._recorder = []
        self._init_data_refs = dict()

    def _comment_(self, text: str):
        self._recorder.append('# %s' % (text, ))

    def _action_(self, action: usage.UsageAction,
                 input_opts: dict, output_opts: dict):
        action_f, action_sig = action.get_action()
        self._update_imports(action_f)

        new_output_opts = {k: k for k in output_opts.keys()}

        t = self._template_action(action_f, input_opts, new_output_opts)
        self._recorder.append(t)

        return new_output_opts

    def _template_action(self, action_f, input_opts, output_opts):
        output_opts = list(output_opts.keys())
        if len(output_opts) == 1:
            output_opts.append('')
        output_vars = ', '.join(output_opts)

        t = '%s = %s(\n' % (output_vars.strip(), action_f.id)
        for k, v in input_opts.items():
            t += '    %s=%s,\n' % (k, v)
        t += ')\n'

        return t

    def _update_imports(self, action_f):
        full_import = action_f.get_import_path()
        import_path, action_api_name = full_import.rsplit('.', 1)
        self._imports.add((import_path, action_api_name))

    def _factory_override(self, ref, factory):
        self._init_data_refs[ref] = factory
        # Don't need to compute anything, so just pass along the ref
        return lambda: ref

    def render(self):
        sorted_imps = sorted(self._imports, key=lambda x: x[0])
        imps = ['from %s import %s\n' % i for i in sorted_imps]
        return '\n'.join(imps + self._recorder)

    def get_example_data(self):
        return {r: f() for r, f in self._init_data_refs.items()}


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
        elif plugin_details[1] == 'pipelines':
            return self._make_spec(name, plugin, ('pipelines',))
        elif plugin_details[1] == 'actions':
            return self._make_spec(name, plugin, ('methods', 'visualizers',
                                                  'pipelines'))
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
            module.pipelines = importlib.import_module('.pipelines',
                                                       package=spec.name)
            module.actions = importlib.import_module('.actions',
                                                     package=spec.name)
        else:
            for action_type in action_types:
                actions = getattr(plugin, action_type)
                for key, value in actions.items():
                    setattr(module, key, value)


sys.meta_path += [QIIMEArtifactAPIImporter()]
