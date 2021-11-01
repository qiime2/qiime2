# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
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


class ArtifactAPIUsageVariable(usage.UsageVariable):
    class quoteless_variable_name:
        def __init__(self, value):
            self.value = value

        def __repr__(self):
            return self.value

    def to_interface_name(self):
        return self.quoteless_variable_name(self.name)


class ArtifactAPIUsage(usage.Usage):
    def __init__(self):
        super().__init__()
        self._reset_state()
        self.global_imports = set()

    def variable_factory(self, name, factory, var_type):
        return ArtifactAPIUsageVariable(
            name,
            factory,
            var_type,
            self,
        )

    def _reset_state(self):
        self.local_imports = set()
        self.recorder = []
        self.init_data_refs = dict()

    def init_artifact(self, name, factory):
        variable = super().init_artifact(name, factory)
        self.init_data_refs[name] = variable
        return variable

    def init_metadata(self, name, factory):
        variable = super().init_metadata(name, factory)
        self.init_data_refs[name] = variable
        return variable

    def merge_metadata(self, name, *variables):
        variable = super().merge_metadata(name, *variables)

        first_md = variables[0].to_interface_name()
        names = [str(r.to_interface_name()) for r in variables[1:]]
        remaining = ', '.join(names)
        t = '%s = %s.merge(%s)\n' % (name, first_md, remaining)
        self.recorder.append(t)

        return variable

    def get_metadata_column(self, name, column_name, variable):
        col_variable = super().get_metadata_column(name, column_name, variable)

        t = '%s = %s.get_column(%r)\n' % (col_variable.to_interface_name(),
                                          variable.to_interface_name(),
                                          column_name)

        self.recorder.append(t)

        return col_variable

    def comment(self, text: str):
        self.recorder.append('# %s' % (text, ))

    def action(self, action, input_opts, output_opts):
        variables = super().action(action, input_opts, output_opts)

        action_f = action.get_action()
        self._update_action_imports(action_f)

        inputs = input_opts.map_variables(lambda v: v.to_interface_name())
        t = self._template_action(action_f, inputs, variables)
        self.recorder.append(t)

        return variables

    def render(self, flush=False):
        sorted_imps = sorted(self.local_imports, key=lambda x: x[0])
        imps = ['from %s import %s\n' % i for i in sorted_imps]
        rendered = '\n'.join(imps + self.recorder)
        if flush:
            self._reset_state()
        return rendered

    def get_example_data(self):
        return {r: f() for r, f in self.init_data_refs.items()}

    def _template_action(self, action_f, input_opts, output_opts):
        outs = [o.to_interface_name() for o in output_opts]
        if len(outs) == 1:
            outs.append('')
        output_vars = ', '.join('%s' % ele for ele in outs)

        t = '%s = %s(\n' % (output_vars.strip(), action_f.id)
        for k, v in input_opts.items():
            if isinstance(v, list):
                t += '    %s=[' % (k,)
                t += ', '.join('%r' % ele for ele in v)
                t += '],\n'
            elif isinstance(v, set):
                t += '    %s={' % (k,)
                t += ', '.join('%r' % ele for ele in sorted(v))
                t += '},\n'
            else:
                t += '    %s=%r,\n' % (k, v)
        t += ')\n'

        return t

    def _update_action_imports(self, action_f):
        full_import = action_f.get_import_path()
        import_path, action_api_name = full_import.rsplit('.', 1)
        import_info = (import_path, action_api_name)
        self._update_imports(import_info)

    def _update_imports(self, import_info):
        if import_info not in self.global_imports:
            self.local_imports.add(import_info)
            self.global_imports.add(import_info)


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
