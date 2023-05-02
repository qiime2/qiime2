from datetime import datetime
import functools
from importlib.metadata import metadata
import pkg_resources
import re
import textwrap
from typing import List, Literal

from q2cli.core.state import get_action_state
from q2cli.core.usage import CLIUsage, CLIUsageVariable
import q2cli.util
from qiime2.core.type import is_semantic_type, is_visualization_type
from qiime2.plugins import ArtifactAPIUsage, ArtifactAPIUsageVariable
from qiime2.sdk.usage import (
    Usage, UsageAction, UsageInputs, UsageOutputNames, UsageOutputs)

from .parse import ProvDAG
from .util import camel_to_snake


class MissingPluginError(Exception):
    """
    an exception class we can use to aggregate missing plugin names
    """
    pass


def _get_action_if_plugin_present(action):
    try:
        return action.get_action()
    except KeyError as e:
        if "No plugin currently registered with id" in (msg := str(e)):
            plugin_id = msg.split()[-1].strip('."\'')
            raise MissingPluginError(
                f"Your QIIME 2 deployment is \n"
                "missing one or more plugins. "
                f"The plugin '{plugin_id}' must be installed to \n"
                "support provenance replay of these Results. "
                "Please install and re-run your command.\n"
                "Many plugins are available at https://library.qiime2.org")


def action_patch(self,
                 action: 'UsageAction',
                 inputs: 'UsageInputs',
                 outputs: 'UsageOutputNames',
                 ) -> 'UsageOutputs':
    """
    A monkeypatch for Usage.action that deals generously with archive versions
    that don't track output-name.

    If there are no output-names to search the signature with, it will attempt
    to search the signature.outputs for a parameter spec with the same
    QIIME type.

    Because our goal in the patched snippet is to assign a usage example type,
    a type-expression match should always return a correct usage example type.
    """
    if not isinstance(action, UsageAction):  # pragma: no cover
        raise ValueError('Invalid value for `action`: expected %r, '
                         'received %r.' % (UsageAction, type(action)))

    if not isinstance(inputs, UsageInputs):  # pragma: no cover
        raise ValueError('Invalid value for `inputs`: expected %r, '
                         'received %r.' % (UsageInputs, type(inputs)))

    if not isinstance(outputs, UsageOutputNames):  # pragma: no cover
        raise ValueError('Invalid value for `outputs`: expected %r, '
                         'received %r.' % (UsageOutputNames,
                                           type(outputs)))

    action_f = _get_action_if_plugin_present(action)

    @functools.lru_cache(maxsize=None)
    def memoized_action():  # pragma: no cover
        execed_inputs = inputs.map_variables(lambda v: v.execute())
        if self.asynchronous:
            return action_f.asynchronous(**execed_inputs).result()
        return action_f(**execed_inputs)

    usage_results = []
    # outputs will be ordered by the `UsageOutputNames` order, not the
    # signature order - this makes it so that the example writer doesn't
    # need to be explicitly aware of the signature order
    for param_name, var_name in outputs.items():
        # param name is not output-name in archive versions without output-name
        try:
            qiime_type = action_f.signature.outputs[param_name].qiime_type
        except KeyError:
            # param_name is often a snake-case qiime2 type, so we can check
            # if the same type still exists in the param spec. If so, use it.
            for (p_name, p_spec) in action_f.signature.outputs.items():
                searchable_type_name = camel_to_snake(str(p_spec.qiime_type))
                if param_name == searchable_type_name:
                    qiime_type = action_f.signature.outputs[p_name].qiime_type
                    break

        if is_visualization_type(qiime_type):
            var_type = 'visualization'
        elif is_semantic_type(qiime_type):
            var_type = 'artifact'
        else:  # pragma: no cover
            raise ValueError('unknown output type: %r' % (qiime_type,))

        def factory(name=param_name):  # pragma: no cover
            results = memoized_action()
            result = getattr(results, name)
            return result

        variable = self._usage_variable(var_name, factory, var_type)
        usage_results.append(variable)

    results = UsageOutputs(outputs.keys(), usage_results)
    cache_info = memoized_action.cache_info
    cache_clear = memoized_action.cache_clear
    # manually graft on cache operations
    object.__setattr__(results, '_cache_info', cache_info)
    object.__setattr__(results, '_cache_reset', cache_clear)
    return results


# NOTE: True monkeypatching happening here. Gross, but the alternative is
# overriding more methods from ReplayCLIUsage and ReplayPythonUsage to point
# to a ReplayUsage subclass
Usage.action = action_patch


def build_header(shebang: str = '', boundary: str = '', copyright: str = '',
                 extra_text: List = []) -> List:
    """
    Writes header copy for all replay outputs, with optional params allowing
    for general utility
    """
    p_lib_md = metadata('qiime2')
    vzn = p_lib_md['Version']
    ts = datetime.now()
    header = []

    if shebang:
        header.append(shebang)

    if boundary:
        header.append(boundary)

    header.extend([
        f"# Auto-generated by provenance_lib v.{vzn} at "
        f"{ts.strftime('%I:%M:%S %p')} on {ts.strftime('%d %b, %Y')}",
        "",
    ])

    if copyright:
        header.extend(copyright)

    header.extend([
        "# For User Support, post to the Community Plugin Support channel of "
        "the QIIME 2",
        "# Forum: https://forum.qiime2.org",
        f"# Documentation/issues: {p_lib_md['Home-page']}",
        "",
        "# UUIDs of all target QIIME 2 Results are shown at the end of "
        "the file",
    ])

    if extra_text:
        header.extend(extra_text)
    if boundary:
        header.append(boundary)
    return header


def build_footer(dag: ProvDAG, boundary: str) -> List:
    """
    Writes footer copy for all outputs
    """
    footer = []
    pairs = []
    uuids = sorted(dag._parsed_artifact_uuids)
    # We can fit two UUIDs on a line, so pair em up
    for idx in range(0, u_len := len(uuids), 2):
        if idx == u_len - 1:
            pairs.append(f"# {uuids[idx]}")
        else:
            pairs.append(f"# {uuids[idx]} \t {uuids[idx + 1]}")

    footer.append(boundary)
    footer.append(
        "# The following QIIME 2 Results were parsed to produce this script:")
    footer.extend(pairs)
    footer.append(boundary)
    footer.append('')
    return footer


class ReplayPythonUsageVariable(ArtifactAPIUsageVariable):
    def to_interface_name(self):
        if self.var_type == 'format':
            return self.name

        parts = {
            'artifact': [self.name],
            'visualization': [self.name, 'viz'],
            'metadata': [self.name, 'md'],
            'column': [self.name],
            # No format here - it shouldn't be possible to make it this far
        }[self.var_type]
        var_name = '_'.join(parts)
        # NOTE: unlike the parent method, this does not guarantee valid python
        # identifiers, because it allows <>. We get more human-readable no-prov
        # node names. Alternately, we could replace < and > with e.g. ___,
        # which is unlikely to occur and is still a valid python identifier
        var_name = re.sub(r'[^a-zA-Z0-9_<>]|^(?=\d)', '_', var_name)
        return self.repr_raw_variable_name(var_name)


class ReplayPythonUsage(ArtifactAPIUsage):
    shebang = '#!/usr/bin/env python'
    header_boundary = '# ' + ('-' * 77)
    copyright = pkg_resources.resource_string(
        __package__, 'assets/copyright_note.txt').decode('utf-8').split('\n')
    how_to = pkg_resources.resource_string(
        __package__, 'assets/python_howto.txt').decode('utf-8').split('\n')

    def __init__(self, enable_assertions: bool = False,
                 action_collection_size: int = 2):
        """Initializer for ReplayPythonUsage
        Identical to ArtifactAPIUsage, but with smaller action_collection_size
        """
        super().__init__()
        self.enable_assertions = enable_assertions
        self.action_collection_size = action_collection_size
        self._reset_state(reset_global_imports=True)

    def _reset_state(self, reset_global_imports=False):
        self.local_imports = set()
        self.header = []
        self.recorder = []
        self.footer = []
        self.init_data_refs = dict()
        if reset_global_imports:
            self.global_imports = set()

    def _template_action(self, action, input_opts, variables):
        """
        Identical to super, but lumps results into `action_results` if
        there are just too many results, and saves results
        """
        action_f = action.get_action()
        if len(variables) > self.action_collection_size or \
                len(action_f.signature.outputs) > 5:
            output_vars = 'action_results'
        else:
            output_vars = self._template_outputs(action, variables)

        plugin_id = action.plugin_id
        action_id = action.action_id
        lines = [
            '%s = %s_actions.%s(' % (output_vars, plugin_id, action_id),
        ]

        all_inputs = (list(action_f.signature.inputs.keys()) +
                      list(action_f.signature.parameters.keys()))
        for k, v in input_opts.items():
            line = ''
            if k not in all_inputs:
                line = self.INDENT + (
                    "# FIXME: The following parameter name was not found in "
                    "your current\n    # QIIME 2 environment. This may occur "
                    "when the plugin version you have\n    # installed does "
                    "not match the version used in the original analysis.\n   "
                    " # Please see the docs and correct the parameter name "
                    "before running.\n")
            line += self._template_input(k, v)
            lines.append(line)

        lines.append(')')

        if len(variables) > self.action_collection_size or \
                len(action.get_action().signature.outputs) > 5:
            for k, v in variables._asdict().items():
                interface_name = v.to_interface_name()
                lines.append('%s = action_results.%s' % (interface_name, k))

        lines.append(
            '# SAVE: comment out the following with \'# \' to skip saving '
            'Results to disk')

        for k, v in variables._asdict().items():
            interface_name = v.to_interface_name()
            lines.append(
                '%s.save(\'%s\')' % (interface_name, interface_name,))

        lines.append('')
        self._add(lines)

    def _template_outputs(self, action, variables):
        """
        Like parent, but allows us to replay an action even when our provenance
        DAG doesn't have a record of all outputs from that action.

        renders outputs we're "not interested in" as _
        """
        output_vars = []
        action_f = action.get_action()

        # need to coax the outputs into the correct order for unpacking
        for output in action_f.signature.outputs:
            try:
                variable = getattr(variables, output)
                output_vars.append(str(variable.to_interface_name()))
            except AttributeError:
                output_vars.append('_')

        if len(output_vars) == 1:
            output_vars.append('')

        return ', '.join(output_vars).strip()

    def init_metadata(self, name, factory, dumped_md_fn: str = ''):
        """
        ArtifactAPIUsage doesn't render Metadata loading, so we do it here.

        dumped_md_fn is an optional parameter used only by
        init_md_from_recorded_md. It allows us to produce scripts that pass
        metadata dumped from provenance to dumped_md_fn
        """
        var = super().init_metadata(name, factory)
        self._update_imports(from_='qiime2', import_='Metadata')
        input_fp = var.to_interface_name()
        if dumped_md_fn:
            lines = [f'{input_fp} = Metadata.load(\'{dumped_md_fn}.tsv\')']
        else:
            self.comment(
                'NOTE: You may substitute already-loaded Metadata for the '
                'following, or cast a pandas.DataFrame to Metadata as needed.'
            )
            lines = [f'{input_fp} = Metadata.load(<your metadata filepath>)']

        self._add(lines)
        return var

    def import_from_format(self, name, semantic_type, variable,
                           view_type=None):
        """
        Identical to super.import_from_format, but writes <your data here>
        instead of import_fp, and saves the result.
        """
        imported_var = Usage.import_from_format(
            self, name, semantic_type, variable, view_type=view_type)

        interface_name = imported_var.to_interface_name()
        import_fp = self.repr_raw_variable_name('<your data here>')

        lines = [
            '%s = Artifact.import_data(' % (interface_name,),
            self.INDENT + '%r,' % (semantic_type,),
            self.INDENT + '%r,' % (import_fp,),
        ]

        if view_type is not None:  # pragma: no cover
            if type(view_type) is not str:
                # Show users where these formats come from when used in the
                # Python API to make things less "magical".
                import_path = super()._canonical_module(view_type)
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

        lines.extend([
            ')',
            '# SAVE: comment out the following with \'# \' to skip saving this'
            ' Result to disk',
            '%s.save(\'%s\')' % (interface_name, interface_name,),
            ''])

        self._update_imports(from_='qiime2', import_='Artifact')
        self._add(lines)

        return imported_var

    def usage_variable(self, name, factory, var_type):
        return ReplayPythonUsageVariable(name, factory, var_type, self)

    class repr_raw_variable_name:
        # allows us to repr col name without enclosing quotes
        # (as in qiime2.qiime2.plugins.ArtifactAPIUsageVariable)
        def __init__(self, value):
            self.value = value

        def __repr__(self):
            return self.value

    def get_metadata_column(self, name, column_name, variable):
        col_variable = Usage.get_metadata_column(
            self, name, column_name, variable)

        to_name = col_variable.to_interface_name()
        from_name = variable.to_interface_name()

        column_name = self.repr_raw_variable_name(column_name)
        lines = ['%s = %s.get_column(%r)' % (to_name, from_name, column_name)]

        self._add(lines)

        return col_variable

    def comment(self, line):
        LINE_LEN = 79
        lines = textwrap.wrap(line, LINE_LEN, break_long_words=False,
                              initial_indent='# ', subsequent_indent='# ')
        lines.append('')
        self._add(lines)

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
        if self.header:
            self.header = self.header + ['']
        if self.footer:
            self.footer = [''] + self.footer
        if sorted_imps:
            sorted_imps = sorted_imps + ['']
        rendered = '\n'.join(
            self.header + sorted_imps + self.recorder + self.footer)
        if flush:
            self._reset_state()
        return rendered

    def build_header(self):
        self.header.extend(
            build_header(self.shebang, self.header_boundary, self.copyright,
                         self.how_to))

    def build_footer(self, dag: ProvDAG):
        self.footer.extend(build_footer(dag, self.header_boundary))


class ReplayCLIUsageVariable(CLIUsageVariable):
    EXT = {
        'artifact': '.qza',
        'visualization': '.qzv',
        'metadata': '',
        'column': '',
        'format': '',
    }

    @property
    def ext(self):
        return self.EXT[self.var_type]

    def to_interface_name(self):
        """
        Like parent, but does not kebab-case metadata. Filepaths are preserved.
        """
        if hasattr(self, '_q2cli_ref'):
            return self._q2cli_ref

        cli_name = '%s%s' % (self.name, self.ext)

        # don't disturb file names, this will break importing where QIIME 2
        # relies on specific filenames being present in a dir
        if self.var_type not in ('format', 'column', 'metadata'):
            cli_name = self.to_cli_name(cli_name)

        return cli_name


class ReplayCLIUsage(CLIUsage):
    shebang = '#!/usr/bin/env bash'
    header_boundary = ('#' * 79)
    copyright = pkg_resources.resource_string(
        __package__, 'assets/copyright_note.txt').decode('utf-8').split('\n')
    set_ex = [
        '# This tells bash to -e exit immediately if a command fails',
        '# and -x show all commands in stdout so you can track progress',
        'set -e -x', '']
    how_to = pkg_resources.resource_string(
        __package__, 'assets/cli_howto.txt').decode('utf-8').split('\n')

    def __init__(self, enable_assertions=False, action_collection_size=None):
        """
        Identical to CLIUsage.init, but creates header & footer attributes
        """
        super().__init__()
        self.header = []
        self.footer = []
        self.enable_assertions = enable_assertions
        self.action_collection_size = action_collection_size

    def usage_variable(self, name, factory, var_type):
        return ReplayCLIUsageVariable(name, factory, var_type, self)

    def _append_action_line(self, signature, param_name, value):
        """
        Like parent but allows replay when recorded parameter names
        are not present in the registered function signatures in the active
        QIIME 2 environment
        """
        param_state = signature.get(param_name)
        if param_state is not None:
            for opt, val in self._make_param(value, param_state):
                line = self.INDENT + opt
                if val is not None:
                    line += ' ' + val
                line += ' \\'
                self.recorder.append(line)
        else:  # no matching param name
            line = self.INDENT + (
                "# FIXME: The following parameter name was not found in "
                "your current\n  # QIIME 2 environment. This may occur "
                "when the plugin version you have\n  # installed does not "
                "match the version used in the original analysis.\n  # "
                "Please see the docs and correct the parameter name "
                "before running.\n")
            cli_name = re.sub('_', '-', param_name)
            line += self.INDENT + '--?-' + cli_name + ' ' + str(value)
            line += ' \\'
            self.recorder.append(line)

    def _make_param(self, value, state):
        """ wrap metadata filenames in <> to force users to replace them """
        if state['metadata'] == 'column':
            value = (f'{value[0]}', *value[1:])
        if state['metadata'] == 'file':
            value = f'{value}'
        return super()._make_param(value, state)

    def import_from_format(self, name, semantic_type, variable,
                           view_type=None):
        """
        Identical to super.import_from_format, but writes --input-path <your
        data here> and follows import block with a blank line
        """
        # We need the super().super() here, so pass self to Usage.import_fr...
        imported_var = Usage.import_from_format(
            self, name, semantic_type, variable, view_type=view_type)

        # in_fp = variable.to_interface_name()
        out_fp = imported_var.to_interface_name()

        lines = [
            'qiime tools import \\',
            self.INDENT + '--type %r \\' % (semantic_type,)
        ]

        if view_type is not None:  # pragma: no cover
            lines.append(
                self.INDENT + '--input-format %s \\' % (view_type,))

        lines += [
            self.INDENT + '--input-path <your data here> \\',
            self.INDENT + '--output-path %s' % (out_fp,),
        ]

        lines.append('')
        self.recorder.extend(lines)

        return imported_var

    def init_metadata(self, name, factory, dumped_md_fn: str = ''):
        """
        Like parent, but appropriately handles filepaths for recorded md fps
        """
        variable = super().init_metadata(name, factory)

        self.init_data.append(variable)

        if dumped_md_fn:
            variable.name = f'"{dumped_md_fn}.tsv"'
        else:
            variable.name = '<your metadata filepath>'

        return variable

    def comment(self, text):
        """
        Identical to parent, but pads comments with an extra newline
        """
        super().comment(text)
        self.recorder.append('')

    def action(self, action, inputs, outputs):
        """
        Overrides parent to fill in missing outputlines from action_f.signature
        Also pads actions with an extra newline
        """
        variables = Usage.action(self, action, inputs, outputs)
        vars_dict = variables._asdict()

        # Get registered collection of output names, so we don't miss any
        # Missing output-names make replay break
        action_f = action.get_action()
        missing_outputs = {}
        for output in action_f.signature.outputs:
            try:
                # If we get a match on output-name, the correct pair is already
                # in vars_dict and we can continue
                getattr(variables, output)
                continue
            except AttributeError:
                # Otherwise, we should add filler values to missing_outputs
                missing_outputs[output] = f'XX_{output}'

        plugin_name = q2cli.util.to_cli_name(action.plugin_id)
        action_name = q2cli.util.to_cli_name(action.action_id)
        self.recorder.append('qiime %s %s \\' % (plugin_name, action_name))

        action_f = action.get_action()
        action_state = get_action_state(action_f)

        ins = inputs.map_variables(lambda v: v.to_interface_name())
        outs = {k: v.to_interface_name() for k, v in vars_dict.items()}
        outs.update(missing_outputs)
        signature = {s['name']: s for s in action_state['signature']}

        for param_name, value in ins.items():
            self._append_action_line(signature, param_name, value)

        max_collection_size = self.action_collection_size
        # NOTE: Losing coverage here until this "lumping" behavior is supported
        if max_collection_size is not None and len(outs) > max_collection_size:
            dir_name = self._build_output_dir_name(plugin_name, action_name)
            self.recorder.append(
                self.INDENT + '--output-dir %s \\' % (dir_name))
            self._rename_outputs(vars_dict, dir_name)
        else:
            for param_name, value in outs.items():
                self._append_action_line(signature, param_name, value)

        self.recorder[-1] = self.recorder[-1][:-2]  # remove trailing \

        self.recorder.append('')
        return variables

    def render(self, flush=False):
        if self.header:
            self.header = self.header + ['']
        if self.footer:
            self.footer = [''] + self.footer
        rendered = '\n'.join(
            self.header + self.set_ex + self.recorder + self.footer)
        if flush:
            self.header = []
            self.footer = []
            self.recorder = []
            self.init_data = []
        return rendered

    def build_header(self):
        self.header.extend(
            build_header(self.shebang, self.header_boundary, self.copyright,
                         self.how_to))

    def build_footer(self, dag: ProvDAG):
        self.footer.extend(build_footer(dag, self.header_boundary))


DRIVER_CHOICES = Literal['python3', 'cli']
SUPPORTED_USAGE_DRIVERS = {
    'python3': ReplayPythonUsage,
    'cli': ReplayCLIUsage,
}
DRIVER_NAMES = list(SUPPORTED_USAGE_DRIVERS.keys())
