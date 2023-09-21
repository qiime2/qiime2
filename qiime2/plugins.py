# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from dataclasses import dataclass
import re
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
    """A specialized implementation for :class:`ArtifactAPIUsage`."""
    # this lets us repr all inputs (including parameters) and have
    # them template out in a consistent manner. without this we would wind
    # up with `foo('my_artifact')` rather than `foo(my_artifact)`.
    class repr_raw_variable_name:
        def __init__(self, value):
            self.value = value

        def __repr__(self):
            return self.value

    def to_interface_name(self):
        if self.var_type == 'format':
            return self.name

        parts = {
            'artifact': [self.name],
            'artifact_collection': [self.name, 'artifact_collection'],
            'visualization': [self.name, 'viz'],
            'visualization_collection': [self.name, 'viz_collection'],
            'metadata': [self.name, 'md'],
            'column': [self.name, 'mdc'],
            # No format here - it shouldn't be possible to make it this far
        }[self.var_type]
        var_name = '_'.join(parts)
        # ensure var_name is a valid python identifier
        var_name = re.sub(r'\W|^(?=\d)', '_', var_name)
        return self.repr_raw_variable_name(var_name)

    def assert_has_line_matching(self, path, expression):
        if not self.use.enable_assertions:
            return

        self.use._update_imports(import_='re')
        name = self.to_interface_name()
        expr = expression

        lines = [
            'hits = sorted(%r._archiver.data_dir.glob(%r))' % (name, path),
            'if len(hits) != 1:',
            self.use.INDENT + 'raise ValueError',
            'target = hits[0].read_text()',
            'match = re.search(%r, target, flags=re.MULTILINE)' % (expr,),
            'if match is None:',
            self.use.INDENT + 'raise AssertionError',
        ]

        self.use._add(lines)

    def assert_output_type(self, semantic_type, key=None):
        if not self.use.enable_assertions:
            return

        name = self.to_interface_name()

        if key:
            name = "%s[%s]" % (name, key)

        lines = [
            'if str(%r.type) != %r:' % (name, str(semantic_type)),
            self.use.INDENT + 'raise AssertionError',
        ]

        self.use._add(lines)


class ArtifactAPIUsage(usage.Usage):
    INDENT = ' ' * 4

    @dataclass(frozen=True)
    class ImporterRecord:
        import_: str
        from_: str = None
        as_: str = None

        def render(self):
            tmpl = 'import %s' % (self.import_,)
            if self.from_ is not None:
                tmpl = 'from %s %s' % (self.from_, tmpl)
            if self.as_ is not None:
                tmpl = '%s as %s' % (tmpl, self.as_)
            return tmpl

    def __init__(self, enable_assertions: bool = False,
                 action_collection_size: int = 3):
        """Constructor for ArtifactAPIUsage

        Warning
        -------
        For SDK use only. Do not use in a written usage example.

        Parameters
        ----------
        enable_assertions : bool
            Whether :class:`qiime2.sdk.usage.UsageVariable` assertions should
            be rendered. Note that these are not executed, rather code that
            would assert is templated by :meth:`render`.
        action_collection_size : int
            The maximum number of outputs to automatically desctructure before
            creating seperate lines with output property access. e.g.
            ``x, y, z = foo()`` vs ``results = foo()`` with ``results.x`` etc.
        """
        super().__init__()
        self.enable_assertions = enable_assertions
        self.action_collection_size = action_collection_size
        self._reset_state(reset_global_imports=True)

    def _reset_state(self, reset_global_imports=False):
        self.local_imports = set()
        self.recorder = []
        self.init_data_refs = dict()
        if reset_global_imports:
            self.global_imports = set()

    def _add(self, lines):
        self.recorder.extend(lines)

    def usage_variable(self, name, factory, var_type):
        return ArtifactAPIUsageVariable(name, factory, var_type, self)

    def render(self, flush: bool = False) -> str:
        """Return a newline-seperated string of Artifact API python code.

        Warning
        -------
        For SDK use only. Do not use in a written usage example.

        Parameters
        ----------
        flush : bool
            Whether to 'flush' the current code. Importantly, this will clear
            the top-line imports for future invocations.
        """
        sorted_imps = sorted(self.local_imports)
        if sorted_imps:
            sorted_imps = sorted_imps + ['']
        rendered = '\n'.join(sorted_imps + self.recorder)
        if flush:
            self._reset_state()
        return rendered

    def init_artifact(self, name, factory):
        variable = super().init_artifact(name, factory)

        var_name = str(variable.to_interface_name())
        self.init_data_refs[var_name] = variable

        return variable

    def init_metadata(self, name, factory):
        variable = super().init_metadata(name, factory)

        var_name = str(variable.to_interface_name())
        self.init_data_refs[var_name] = variable

        return variable

    def init_artifact_collection(self, name, factory):
        variable = super().init_artifact_collection(name, factory)

        var_name = str(variable.to_interface_name())
        self.init_data_refs[var_name] = variable

        return variable

    def init_format(self, name, factory, ext=None):
        if ext is not None:
            name = '%s.%s' % (name, ext.lstrip('.'))

        variable = super().init_format(name, factory, ext=ext)

        var_name = str(variable.to_interface_name())
        self.init_data_refs[var_name] = variable

        return variable

    def import_from_format(self, name, semantic_type,
                           variable, view_type=None):
        imported_var = super().import_from_format(
            name, semantic_type, variable, view_type=view_type)

        interface_name = imported_var.to_interface_name()
        import_fp = variable.to_interface_name()

        lines = [
            '%s = Artifact.import_data(' % (interface_name,),
            self.INDENT + '%r,' % (semantic_type,),
            self.INDENT + '%r,' % (import_fp,),
        ]

        if view_type is not None:
            if type(view_type) is not str:
                # Show users where these formats come from when used in the
                # Python API to make things less "magical".
                import_path = _canonical_module(view_type)
                view_type = view_type.__name__
                if import_path is not None:
                    self._update_imports(from_=import_path,
                                         import_=view_type)
                else:
                    # May be in scope already, but something is quite wrong at
                    # this point, so assume the plugin_manager is sufficiently
                    # informed.
                    view_type = repr(view_type)
            else:
                view_type = repr(view_type)

            lines.append(self.INDENT + '%s,' % (view_type,))

        lines.append(')')

        self._update_imports(from_='qiime2', import_='Artifact')
        self._add(lines)

        return imported_var

    def merge_metadata(self, name, *variables):
        variable = super().merge_metadata(name, *variables)

        first_var, remaining_vars = variables[0], variables[1:]
        first_md = first_var.to_interface_name()

        names = [str(r.to_interface_name()) for r in remaining_vars]
        remaining = ', '.join(names)
        var_name = variable.to_interface_name()

        lines = ['%r = %r.merge(%s)' % (var_name, first_md, remaining)]

        self._add(lines)

        return variable

    def get_metadata_column(self, name, column_name, variable):
        col_variable = super().get_metadata_column(name, column_name, variable)

        to_name = col_variable.to_interface_name()
        from_name = variable.to_interface_name()

        lines = ['%s = %s.get_column(%r)' % (to_name, from_name, column_name)]

        self._add(lines)

        return col_variable

    def view_as_metadata(self, name, from_variable):
        to_variable = super().view_as_metadata(name, from_variable)

        from_name = from_variable.to_interface_name()
        to_name = to_variable.to_interface_name()

        lines = ['%r = %r.view(Metadata)' % (to_name, from_name)]

        self._update_imports(from_='qiime2', import_='Metadata')
        self._add(lines)

        return to_variable

    def peek(self, variable):
        var_name = variable.to_interface_name()

        lines = []
        for attr in ('uuid', 'type', 'format'):
            lines.append('print(%r.%s)' % (var_name, attr))

        self._add(lines)

        return variable

    def comment(self, text):
        lines = ['# %s' % (text,)]

        self._add(lines)

    def help(self, action):
        action_name = self._plugin_import_as_name(action)

        # TODO: this isn't pretty, but it gets the job done
        lines = ['help(%s.%s.__call__)' % (action_name, action.action_id)]

        self._add(lines)

    def action(self, action, input_opts, output_opts):
        variables = super().action(action, input_opts, output_opts)

        self._plugin_import_as_name(action)

        inputs = input_opts.map_variables(lambda v: v.to_interface_name())
        self._template_action(action, inputs, variables)

        return variables

    def get_example_data(self):
        return {r: v.execute() for r, v in self.init_data_refs.items()}

    def _plugin_import_as_name(self, action):
        action_f = action.get_action()
        full_import = action_f.get_import_path()
        base, _, _ = full_import.rsplit('.', 2)
        as_ = '%s_actions' % (action.plugin_id,)
        self._update_imports(import_='%s.actions' % (base,), as_=as_)
        return as_

    def _template_action(self, action, input_opts, variables):
        if len(variables) > self.action_collection_size:
            output_vars = 'action_results'
        else:
            output_vars = self._template_outputs(action, variables)

        plugin_id = action.plugin_id
        action_id = action.action_id
        lines = [
            '%s = %s_actions.%s(' % (output_vars, plugin_id, action_id),
        ]

        for k, v in input_opts.items():
            line = self._template_input(k, v)
            lines.append(line)

        lines.append(')')

        if len(variables) > self.action_collection_size:
            for k, v in variables._asdict().items():
                var_name = v.to_interface_name()
                lines.append('%s = action_results.%s' % (var_name, k))

        self._add(lines)

    def _template_outputs(self, action, variables):
        output_vars = []
        action_f = action.get_action()

        # need to coax the outputs into the correct order for unpacking
        for output in action_f.signature.outputs:
            variable = getattr(variables, output)
            output_vars.append(str(variable.to_interface_name()))

        if len(output_vars) == 1:
            output_vars.append('')

        return ', '.join(output_vars).strip()

    def _template_input(self, input_name, value):
        if isinstance(value, list):
            t = ', '.join(repr(el) for el in value)
            return self.INDENT + '%s=[%s],' % (input_name, t)

        if isinstance(value, set):
            t = ', '.join(repr(el) for el in sorted(value, key=str))
            return self.INDENT + '%s={%s},' % (input_name, t)

        return self.INDENT + '%s=%r,' % (input_name, value)

    def _update_imports(self, import_, from_=None, as_=None):
        import_record = self.ImporterRecord(
            import_=import_, from_=from_, as_=as_)

        if as_ is not None:
            self.namespace.add(as_)
        else:
            self.namespace.add(import_)

        rendered = import_record.render()
        if rendered not in self.global_imports:
            self.local_imports.add(rendered)
            self.global_imports.add(rendered)


def _canonical_module(obj):
    last_module = None
    module_str = obj.__module__
    parts = module_str.split('.')
    while parts:
        try:
            module = importlib.import_module('.'.join(parts))
        except ModuleNotFoundError:
            return last_module
        if not hasattr(module, obj.__name__):
            return last_module

        last_module = '.'.join(parts)
        parts.pop()

    return None


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
