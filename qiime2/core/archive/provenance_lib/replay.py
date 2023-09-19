import bibtexparser as bp
from bibtexparser.bwriter import BibTexWriter
import networkx as nx
import os
import pathlib
import pkg_resources
import shutil
import tempfile
from collections import UserDict
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Set, Union

from .archive_parser import ProvNode
from .parse import ProvDAG
from .usage_drivers import (
    DRIVER_CHOICES, DRIVER_NAMES, SUPPORTED_USAGE_DRIVERS, Usage,
    build_header, build_footer
)
from ..provenance import MetadataInfo

from qiime2.sdk import PluginManager
from qiime2.sdk.usage import UsageVariable
from qiime2.sdk.util import camel_to_snake


@dataclass
class ReplayConfig():
    '''
    Dataclass that stores various user-selected configuration options and
    other bits of information relevant to provenance replay.

    Parameters
    ----------
    use : Usage
        The usage driver to be used for provenance replay.
    dump_recorded_metadata : bool
        If True, replay should write the metadata recorded in provenance to
        disk in .tsv format.
    use_recorded_metadata : bool
        If True, replay should use the metadata recorded in provenance.
    pm : PluginManager
        The active instance of the QIIME 2 PluginManager.
    md_context_has_been_printed : bool
        A flag set by default and used internally, allows context to be
        printed once and only once.
    no_provenance_context_has_been_printed : bool
        Indicates whether the no-provenance context documentation has been
        printed.
    header : bool
        If True, an introductory how-to header should be rendered in the
        script.
    verbose : bool
        If True, progress is reported to stdout.
    md_out_fp : str
        The directory where caputred metadata should be written.
    '''
    use: Usage
    dump_recorded_metadata: bool = True
    use_recorded_metadata: bool = False
    pm: PluginManager = PluginManager()
    md_context_has_been_printed: bool = False
    no_provenance_context_has_been_printed: bool = False
    header: bool = True
    verbose: bool = False
    md_out_fp: str = ''


@dataclass
class ActionCollections():
    '''
    std_actions are all normal provenance-tracked q2 actions, arranged as:

    {
        <action_id>: {
            <output_node_uuid>: 'output_name',
            <output_node_2_uuid:> 'output_name_2'
        },
        <action_2_id> : ...
     }

    no_provenance_nodes can't be organized by action, and in some cases we
    don't know anything but UUID for them, so we can fit these in a list.
    '''
    std_actions: Dict[str, Dict[str, str]] = field(default_factory=dict)
    no_provenance_nodes: List[str] = field(default_factory=list)


class UsageVarsDict(UserDict):
    '''
    Mapping of uuid -> unique variable name.

    A dict where values are also unique. Used here as a UUID-queryable
    "namespace" of strings that can be passed to usage drivers for rendering
    into unique variable names. Non-unique values cause namespace collisions.

    For consistency and simplicity, all str values are suffixed with _n when
    added, such that n is some int. When potentially colliding values are
    added, n is incremented as needed until collision is avoided.
    UsageVarsDicts mutate ALL str values they receive.

    Best practice is generally to add the UUID: variable-name pair to this,
    create the usage variable using the stored name,
    then store the usage variable in a separate {UUID: UsageVar}. This
    ensures that UsageVariable.name is unique, preventing namespace collisions.
    NamespaceCollections (below) exist to group these related structures.

    Note: it's not necessary (and may break the mechanism of uniqueness here)
    to maintain parity between variable names in this namespace and in the
    usage variable store. The keys in both stores, however, must match.
    '''
    def __setitem__(self, key: str, item: str) -> None:
        unique_item = self._uniquify(item)
        return super().__setitem__(key, unique_item)

    def _uniquify(self, var_name: str) -> str:
        '''
        Appends _<some int> to var_name, such that the returned name won't
        collide with any variable-name values that already exist in the dict.

        Parameters
        ----------
        var_name : str
            The variable name to make unique.

        Returns
        -------
        str
            The unique integer-appended variable name.
        '''
        some_int = 0
        unique_name = f'{var_name}_{some_int}'
        values = self.data.values()

        # no-prov nodes are stored with angle brackets around them, but
        # those brackets shouldn't be considered on uniqueness check
        while unique_name in values or f'<{unique_name}>' in values:
            some_int += 1
            unique_name = f"{var_name}_{some_int}"
        return unique_name

    def get_key(self, value: str):
        '''
        Given some value in the dict, returns its key.

        Parameters
        ----------
        value : str
            The value to query.

        Returns
        -------
        str
            The key (uuid) corresponding to the value (variable name).

        Raises
        ------
        KeyError
            If the value is not found.
        '''
        for key, val in self.items():
            if value == val:
                return key
        raise KeyError(f'passed value \'{value}\' does not exist in the dict.')

    def wrap_val_in_angle_brackets(self, key: str):
        '''
        Wraps the variable name pointed to by `key` in brackets `<>`.

        Parameters
        ----------
        key : str
            The key (uuid) of the variable name to wrap in brackets.

        Notes
        -----
        If this is accidentally run twice, it will break things.
        '''
        super().__setitem__(key, f'<{self.data[key]}>')


@dataclass
class NamespaceCollections:
    '''
    A dataclass collection of objects that each track some useful bit of
    information relevant to usage namespaces.

    Attributes
    ----------
    usg_var_namespace : UsageVarsDict
        A uuid -> variable-name mapping that ensures that variable names remain
        uniue in the dictionary.
    usg_vars : dict
        A uuid -> UsageVariable mapping. The names of the UsageVariables here
        do not necessarily match those in `usg_var_namespace`.
    action_namespace : set of str
        A collection of unique action strings that look like
        `{plugin}_{action}_{sequential int}`.
    '''
    usg_var_namespace: UsageVarsDict = field(default_factory=UsageVarsDict)
    usg_vars: Dict[str, UsageVariable] = field(default_factory=dict)
    action_namespace: Set[str] = field(default_factory=set)


def replay_provenance(
    payload: Union[str, ProvDAG],
    out_fp: str,
    usage_driver: DRIVER_CHOICES = 'python3',
    validate_checksums: bool = True,
    parse_metadata: bool = True,
    recurse: bool = False,
    use_recorded_metadata: bool = False,
    suppress_header: bool = False,
    verbose: bool = False,
    dump_recorded_metadata: bool = True,
    md_out_fp: str = ''
):
    '''
    Renders usage examples describing a ProvDAG, producing an interface-
    specific executable.

    ProvDAG inputs retain their original config values. The
    `validate_checksums`, `parse_metadata`, `recurse`, and `verbose` parameters
    are disregarded if the payload is a ProvDAG.

    Parameters
    ----------
    payload : str or ProvDAG
        A filepath to an artifact or directory containing artifacts, or the
        ProvDAG to be parsed.
    out_fp : str
        The filepath at which to write the rendered executable.
    usage_driver : one of 'python3', 'cli'
        The interface to target.
    validate_checksums : bool
        Whether to perform checksum validation on the input artifact.
    parse_metadata : bool
        Whether to parse study metadata recorded in provenance.
    recurse : bool
        Whether to recursively parse nested directories containing artifacts.
    use_recorded_metadata : bool
        Whether to use the metadata recorded in provenance.
    suppress_header : bool
        Whether to forgo rendering the header and footer that are included
        by default in replay scripts.
    verbose : bool
        Whether to print status messages during processing.
    dump_recorded_metadata : bool
        Whether to write the metadata recorded in provenance to disk.
    md_out_fp : str
        The directory in which to write the recorded metadata if desired.
    '''
    if type(payload) is ProvDAG:
        parse_metadata = payload.cfg.parse_study_metadata

    if not parse_metadata:
        if use_recorded_metadata:
            raise ValueError(
                'Metadata not parsed for replay. Re-run with parse_metadata, '
                'or set use_recorded_metadata to False.'
            )
        if dump_recorded_metadata:
            raise ValueError(
                'Metadata not parsed, so cannot be written to disk. Re-run '
                'with parse_metadata, or set dump_recorded_metadata to False.'
            )
        if md_out_fp:
            raise ValueError(
                'Metadata not parsed, so cannot be written to disk. Re-run '
                'with parse_metadata, or do not pass a metadata output '
                'filepath argument.'
            )

    if use_recorded_metadata and not dump_recorded_metadata:
        raise NotImplementedError(
            'In order to produce a replay script that uses metadata '
            'captured in provenance, that metadata must first be written to '
            'disk. Re-run with dump-recorded-metadata set to True, or '
            'use-recorded-metadata set to False.'
        )

    dag = ProvDAG(
        payload, validate_checksums, parse_metadata, recurse, verbose
    )
    cfg = ReplayConfig(
        use=SUPPORTED_USAGE_DRIVERS[usage_driver](),
        use_recorded_metadata=use_recorded_metadata,
        dump_recorded_metadata=dump_recorded_metadata,
        verbose=verbose, md_out_fp=md_out_fp
    )

    build_usage_examples(dag, cfg)
    if not suppress_header:
        cfg.use.build_header()
        cfg.use.build_footer(dag)

    if cfg.dump_recorded_metadata:
        print('metadata written to recorded_metadata/')

    output = cfg.use.render(flush=True)
    with open(out_fp, mode='w') as out_fh:
        out_fh.write(output)


def group_by_action(dag: ProvDAG, nodes: Iterator[str]) -> ActionCollections:
    '''
    This groups the nodes from a DAG by action, returning an ActionCollections
    aggregating the outputs related to each action.

    Takes an iterator of UUIDs, allowing us to influence the ordering of the
    grouping.

    In cases where a captured output_name is unavailable, we substitute the
    output data's Semantic Type, snake-cased because it will be used as
    a variable name if this data is rendered by ArtifactAPIUsage.

    Parameters
    ----------
    dag : ProvDAG
        The dag representation of parsed provenance.
    nodes : iterator of str
        An iterator over node uuids.

    Returns
    -------
    ActionCollections
        The outputs grouped by action.
    '''
    actions = ActionCollections()
    for node_id in nodes:
        if dag.node_has_provenance(node_id):
            node = dag.get_node_data(node_id)
            action_id = node.action._execution_details['uuid']

            output_name = node.action.output_name
            if output_name is None:
                output_name = camel_to_snake(node.type)

            try:
                actions.std_actions[action_id].update({node_id: output_name})
            except KeyError:
                actions.std_actions[action_id] = {node_id: output_name}
        else:
            actions.no_provenance_nodes.append(node_id)

    return actions


def build_usage_examples(dag: ProvDAG, cfg: ReplayConfig):
    '''
    Builds a chained usage example representing the analysis `dag`.

    Parameters
    ----------
    dag : ProvDAG
       The dag representation of parsed provenance.
    cfg : ReplayConfig
        Replay configuration options.
    '''
    usg_ns = NamespaceCollections()
    sorted_nodes = nx.topological_sort(dag.collapsed_view)
    actions = group_by_action(dag, sorted_nodes)

    for node_id in actions.no_provenance_nodes:
        node = dag.get_node_data(node_id)
        build_no_provenance_node_usage(node, node_id, usg_ns, cfg)

    for action_id in (std_actions := actions.std_actions):
        # we are replaying actions not nodes, so any associated node works
        some_node_id = next(iter(std_actions[action_id]))
        node = dag.get_node_data(some_node_id)
        if node.action.action_type == 'import':
            build_import_usage(node, usg_ns, cfg)
        else:
            build_action_usage(node, usg_ns, std_actions, action_id, cfg)


def build_no_provenance_node_usage(
    node: Optional[ProvNode],
    uuid: str,
    ns: NamespaceCollections,
    cfg: ReplayConfig
):
    '''
    Given a ProvNode (with no provenance), make sure comments will be rendered
    explaining this, add an empty usage variable to the namespace and log the
    node. Returns nothing, modifying the passed usage instance in place.

    Parameters
    ----------
    node : ProvNode or None
        Either a no-provenance node, or None indicating that only `uuid`
        is available.
    uuid : str
        The uuid of the node/result.
    ns : NamespaceCollections
        Info tracking usage namespaces.
    cfg : ReplayConfig
        Replay configuration options. Contains the modified usage driver.
    '''
    if not cfg.no_provenance_context_has_been_printed:
        cfg.no_provenance_context_has_been_printed = True
        cfg.use.comment(
            'One or more nodes have no provenance, so full replay is '
            'impossible. Any commands we were able to reconstruct have been '
            'rendered, with the string descriptions below replacing actual '
            'inputs.'
        )
        cfg.use.comment(
            'Original Node ID                       String Description'
        )
    if node is None:
        # the node is a !no-provenance input and we have only UUID
        var_name = 'no-provenance-node'
    else:
        var_name = camel_to_snake(node.type)
    ns.usg_var_namespace.update({uuid: var_name})
    ns.usg_var_namespace.wrap_val_in_angle_brackets(uuid)

    # make a usage variable for downstream consumption
    empty_var = cfg.use.usage_variable(
        ns.usg_var_namespace[uuid], lambda: None, 'artifact'
    )
    ns.usg_vars.update({uuid: empty_var})

    # log the no-prov node
    cfg.use.comment(f"{uuid}   {ns.usg_vars[uuid].to_interface_name()}")


def build_import_usage(
        node: ProvNode, ns: NamespaceCollections, cfg: ReplayConfig
):
    '''
    Given a ProvNode, adds an import usage example for it, roughly
    resembling the below. Returns nothing, modifying the passed usage instance
    in place.

    raw_seqs = use.init_format('raw_seqs', lambda: None, ext='fastq.gz')
    imported_seqs = use.import_from_format(
        'emp_single_end_sequences',
        'EMPSingleEndSequences',
        raw_seqs
    )
    The `lambda: None` is a placeholder for some actual data factory,
    and should not impact the rendered usage.

    Parameters
    ----------
    node : ProvNode
        The imported node of interest.
    ns : NamespaceCollections
        Info tracking usage namespaces.
    cfg : ReplayConfig
        Replay configuration options. Contains the modified usage driver.
    '''
    format_id = node._uuid + '_f'
    ns.usg_var_namespace.update({format_id: camel_to_snake(node.type) + '_f'})
    format_for_import = cfg.use.init_format(
        ns.usg_var_namespace[format_id], lambda: None
    )

    ns.usg_var_namespace.update({node._uuid: camel_to_snake(node.type)})
    use_var = cfg.use.import_from_format(
        ns.usg_var_namespace[node._uuid], node.type, format_for_import
    )
    ns.usg_vars.update({node._uuid: use_var})


def build_action_usage(
    node: ProvNode,
    ns: NamespaceCollections,
    std_actions: Dict[str, Dict[str, str]],
    action_id: str,
    cfg: ReplayConfig
):
    '''
    Adds an action usage example to `use` for some ProvNode.
    Returns nothing, modifying the passed usage instance in place.

    use.action(
        use.UsageAction(plugin_id='diversity_lib',
                        action_id='pielou_evenness'),
        use.UsageInputs(table=ft),
        use.UsageOutputNames(vector='pielou_vector')
    )

    Parameters
    ----------
    node : ProvNode
        The node the creating action of which is of interest.
    ns : NamespaceCollections
        Info tracking usage namespaces.
    std_actions : dict
        Expalained in ActionCollections.
    action_id : str
        The uuid of the action.
    cfg : ReplayConfig
        Replay configuration options. Contains the modified usage driver.
    '''
    command_specific_md_context_has_been_printed = False
    plugin = node.action.plugin
    action = node.action.action_name
    plg_action_name = uniquify_action_name(plugin, action, ns.action_namespace)

    inputs = _collect_action_inputs(ns, node)

    # Process outputs before params so we can access the unique output name
    # from the namespace when dumping metadata to files below
    raw_outputs = std_actions[action_id].items()
    outputs = _uniquify_output_names(ns, raw_outputs)

    for param_name, param_val in node.action.parameters.items():
        # We can currently assume that None arguments are only passed to params
        # as default values, so we can skip these parameters entirely in replay
        if param_val is None:
            continue

        if isinstance(param_val, MetadataInfo):
            unique_md_id = ns.usg_var_namespace[node._uuid] + '_' + param_name
            ns.usg_var_namespace.update(
                {unique_md_id: camel_to_snake(param_name)}
            )
            md_fn = ns.usg_var_namespace[unique_md_id]
            if cfg.dump_recorded_metadata:
                md_with_ext = md_fn + '.tsv'
                dump_recorded_md_file(
                    cfg, node, plg_action_name, param_name, md_with_ext
                )

            if cfg.use_recorded_metadata:
                # the local dir and fp where md will be saved (if at all) is:
                md_fn = f'{plg_action_name}/{md_fn}'
                md = init_md_from_recorded_md(
                    node,
                    param_name,
                    unique_md_id,
                    ns.usg_var_namespace,
                    cfg,
                    md_fn
                )
            else:
                if not cfg.md_context_has_been_printed:
                    cfg.md_context_has_been_printed = True
                    cfg.use.comment(
                        "Replay attempts to represent metadata inputs "
                        "accurately, but metadata .tsv files are merged "
                        "automatically by some interfaces, rendering "
                        "distinctions between file inputs invisible in "
                        "provenance. We output the recorded metadata to disk "
                        "to enable visual inspection.")

                if not command_specific_md_context_has_been_printed:
                    if cfg.md_out_fp:
                        fp = f'{cfg.md_out_fp}/{plg_action_name}'
                    else:
                        fp = f'./recorded_metadata/{plg_action_name}/'

                    cfg.use.comment(
                        "The following command may have received additional "
                        "metadata .tsv files. To confirm you have covered "
                        "your metadata needs adequately, review the original "
                        f"metadata, saved at '{fp}'")

                if not param_val.input_artifact_uuids:
                    md = init_md_from_md_file(
                        node,
                        param_name,
                        unique_md_id,
                        ns.usg_var_namespace,
                        cfg
                    )
                else:
                    md = init_md_from_artifacts(param_val, ns, cfg)

            param_val = md

        inputs.update({param_name: param_val})

    usg_var = cfg.use.action(
        cfg.use.UsageAction(plugin_id=plugin, action_id=action),
        cfg.use.UsageInputs(**inputs),
        cfg.use.UsageOutputNames(**outputs)
    )

    # write the usage vars into the UsageVars dict so we can use em downstream
    for res in usg_var:
        uuid_key = ns.usg_var_namespace.get_key(value=res.name)
        ns.usg_vars[uuid_key] = res


def _collect_action_inputs(ns: NamespaceCollections, node: ProvNode) -> dict:
    '''
    Returns a dict containing the action Inputs for a ProvNode.
    Dict structure: {input_name: input_var} or {input_name: [input_var1, ...]}.

    Parameters
    ----------
    ns : NamespaceCollections
        Info tracking usage namespaces.
    node : ProvNode
        The node the creating action of which's inputs are of interest.

    Returns
    -------
    dict
        Mapping input names to their corresponding usage variables.
    '''
    inputs_dict = {}
    for input_name, uuids in node.action.inputs.items():
        # Some optional inputs take None as a default
        if uuids is not None:
            if type(uuids) is str:
                uuid = uuids
                inputs_dict.update({input_name: ns.usg_vars[uuid]})
            else:  # it's a collection
                input_vars = []
                for uuid in uuids:
                    input_vars.append(ns.usg_vars[uuid])
                inputs_dict.update({input_name: input_vars})
    return inputs_dict


def _uniquify_output_names(
    ns: NamespaceCollections, raw_outputs: dict
) -> dict:
    '''
    Returns a dict containing the uniquified output names from a ProvNode.
    Dict structure: {output_name: uniquified_output_name}.

    Parameters
    ----------
    ns : NamespaceCollections
        Info tracking usage namespaces.
    raw_outputs : dict
        Mapping of node uuid to output-name as seen in action.yaml.

    Returns
    -------
    dict
        Mapping of original output-name to output-name after being made unique.
    '''
    outputs = {}
    for uuid, output_name in raw_outputs:
        ns.usg_var_namespace.update({uuid: output_name})
        uniquified_output_name = ns.usg_var_namespace[uuid]
        outputs.update({output_name: uniquified_output_name})
    return outputs


def init_md_from_recorded_md(
    node: ProvNode,
    param_name: str,
    md_id: str,
    ns: UsageVarsDict,
    cfg: ReplayConfig,
    md_fn: str
) -> UsageVariable:
    '''
    Initializes and returns a Metadata UsageVariable from Metadata parsed
    from provenance.

    Parameters
    ----------
    node : ProvNode
        The node the creating action of which was passed metadata.
    param_name : str
        The name of the parameter to which metadata was passed.
    md_id : str
        Looks like: f'{node uuid}_{param name}'.
    ns : UsageVarsDict
        Mapping of uuid -> unique variable name.
    cfg : ReplayConfig
        Replay configuration options. Contains the executing usage driver.
    md_fn : str
        Looks like: f'{plugin}_{action}_{counter}/{unique param name}.

    Returns
    -------
    UsageVariable
        Of type metadata or metadata column.

    Raises
    ------
    ValueError
        If the node has no metadata.
    '''
    if not node.metadata:
        raise ValueError(
            'This function should only be called if the node has metadata.'
        )

    md_df = node.metadata[param_name]

    def factory():
        from qiime2 import Metadata
        return Metadata(md_df)

    cwd = pathlib.Path.cwd()
    if cfg.md_out_fp:
        fn = str(cwd / cfg.md_out_fp / md_fn)
    else:
        fn = str(cwd / 'recorded_metadata' / md_fn)

    md = cfg.use.init_metadata(ns[md_id], factory, dumped_md_fn=fn)
    plugin = node.action.plugin
    action = node.action.action_name
    if param_is_metadata_column(cfg, param_name, plugin, action):
        mdc_id = node._uuid + '_mdc'
        mdc_name = ns[md_id] + '_mdc'
        ns.update({mdc_id: mdc_name})
        md = cfg.use.get_metadata_column(ns[mdc_id], '<column name>', md)
    return md


def init_md_from_md_file(
    node: ProvNode,
    param_name: str,
    md_id: str,
    ns: UsageVarsDict,
    cfg: ReplayConfig
) -> UsageVariable:
    '''
    Initializes and returns a Metadata UsageVariable with no real data,
    mimicking a user passing md as a .tsv file.

    Parameters
    ----------
    node : ProvNode
        The node the creating action of which was passed a metadata file.
    param_name : str
        The parameter name to which the metadata file was passed.
    md_id : str
        Looks like: f'{node uuid}_{param name}'.
    ns : UsageVarsDict
        Mapping of uuid -> unique variable name.
    cfg : ReplayConfig
        Replay configuration options. Contains the executing usage driver.

    Returns
    -------
    UsageVariable
        Of type metadata or metadata column.
    '''
    plugin = node.action.plugin
    action = node.action.action_name
    md = cfg.use.init_metadata(ns[md_id], lambda: None)
    if param_is_metadata_column(cfg, param_name, plugin, action):
        mdc_id = node._uuid + '_mdc'
        mdc_name = ns[md_id] + '_mdc'
        ns.update({mdc_id: mdc_name})
        md = cfg.use.get_metadata_column(ns[mdc_id], '<column name>', md)
    return md


def init_md_from_artifacts(
        md_inf: MetadataInfo, ns: NamespaceCollections, cfg: ReplayConfig
) -> UsageVariable:
    '''
    Initializes and returns a Metadata UsageVariable with no real data,
    mimicking a user passing one or more QIIME 2 Artifacts as metadata.

    We expect these usage vars are already in the namespace as artifacts if
    we're reading them in as metadata.

    Parameters
    ----------
    md_inf : MetadataInfo
        Named tuple with fields `input_artifact_uuids` which is a list of
        uuids and `relative_fp` which is the filename of the metadata file.
        These are parsed from a !metadata tag in action.yaml.
    ns : NamespaceCollections
        Info tracking usage namespaces.
    cfg: ReplayConfig
        Replay configuration options. Contains the executing usage driver.

    Returns
    -------
    UsageVariable
        Of type metadata.

    Raises
    ------
    ValueError
        If no input artifact uuids are present in MetadataInfo.
    '''
    if not md_inf.input_artifact_uuids:
        raise ValueError(
            'This funtion should not be used if '
            'MetadataInfo.input_artifact_uuids is empty.'
        )

    md_files_in = []
    for artifact_uuid in md_inf.input_artifact_uuids:
        amd_id = artifact_uuid + '_a'
        var_name = ns.usg_vars[artifact_uuid].name + '_a'
        if amd_id not in ns.usg_var_namespace:
            ns.usg_var_namespace.update({amd_id: var_name})
            art_as_md = cfg.use.view_as_metadata(
                ns.usg_var_namespace[amd_id], ns.usg_vars[artifact_uuid]
            )
            ns.usg_vars.update({amd_id: art_as_md})
        else:
            art_as_md = ns.usg_vars[amd_id]
        md_files_in.append(art_as_md)

    if len(md_inf.input_artifact_uuids) > 1:
        # we can't uniquify this normally, because one uuid can be merged with
        # combinations of others
        merge_id = '-'.join(md_inf.input_artifact_uuids)
        ns.usg_var_namespace.update({merge_id: 'merged_artifacts'})
        merged_md = cfg.use.merge_metadata(
            ns.usg_var_namespace[merge_id], *md_files_in
        )
        ns.usg_vars.update({merge_id: merged_md})
    return art_as_md


def dump_recorded_md_file(
    cfg: ReplayConfig,
    node: ProvNode,
    action_name: str,
    md_id: str,
    fn: str
):
    '''
    Writes one metadata DataFrame pointed to by `md_id` to a .tsv file.
    Each action gets its own directory containing relevant md files.

    Raises a ValueError if the node has no metadata

    Parameters
    ----------
    cfg : ReplayConfig
        Replay configuration options. Contains the executing usage driver.
    node : ProvNode
        The node the creating action of which recorded metadata as input. Used
        here only to ensure that metadata was in fact recorded.
    action_name : str
        Looks like: f'{plugin}_{action}_{counter}'.
    md_id : str
        Looks like: f'{node uuid}_{param name}'.
    fn : str
        Looks like: f'{unique param name}.tsv'.

    Raises
    ------
    ValueError
        If the passed node does not have metadata in its creating action.
    '''
    if node.metadata is None:
        raise ValueError(
            'This function should only be called if the node has metadata.'
        )

    if cfg.md_out_fp:
        md_out_fp_base = pathlib.Path(cfg.md_out_fp)
    else:
        cwd = pathlib.Path.cwd()
        md_out_fp_base = cwd / 'recorded_metadata'

    action_dir = md_out_fp_base / action_name
    action_dir.mkdir(parents=True, exist_ok=True)

    md_df = node.metadata[md_id]
    out_fp = action_dir / (fn)
    md_df.to_csv(out_fp, sep='\t', index=False)


def param_is_metadata_column(
    cfg: ReplayConfig, param: str, plugin: str, action: str
) -> bool:
    '''
    Returns True if the parameter name `param` is registered as a
    MetadataColumn.

    Parameters
    ----------
    cfg : ReplayConfig
        Replay configuration options. Contains the plugin manager object.
    param : str
        The name of the parameter of interest.
    plugin : str
        The plugin that the relevant action belongs to.
    action : str
        The action that has the parameter of interest.

    Returns
    -------
    bool
        Indicating whether the parameter of interest is a MetadataColumn.

    Raises
    ------
    KeyError
        - If the plugin of interest is not registered with the plugin manager.
        - If the action of interest is not registered with the plugin.
        - If the parameter is not in the signature of the action.
    '''
    plugin = cfg.pm.get_plugin(id=plugin)

    try:
        action_f = plugin.actions[action]
    except KeyError:
        raise KeyError(
            f'No action registered with name {action} in plugin {plugin}.'
        )

    try:
        param_spec = action_f.signature.parameters[param]
    except KeyError:
        raise KeyError(
            f'No parameter registered with name {param} in action {action}.'
        )

    # HACK, but it works without relying on Q2's type system
    return 'MetadataColumn' in str(param_spec.qiime_type)


def uniquify_action_name(
    plugin: str, action: str, action_namespace: Set[str]
) -> str:
    '''
    Creates a unique name by concatenating plugin, action, and a counter,
    and adds this name to action_ns before returning it.

    Parameters
    ----------
    plugin : str
        The name of the plugin.
    action : str
        The name of the action.
    action_namespace : set of str
        The collection of unqiue action names.

    Returns
    -------
    str
        The unique action name.
    '''
    counter = 0
    plg_action_name = f'{plugin}_{action}_{counter}'
    while plg_action_name in action_namespace:
        counter += 1
        plg_action_name = f'{plugin}_{action}_{counter}'
    action_namespace.add(plg_action_name)
    return plg_action_name


def collect_citations(
    dag: ProvDAG, deduplicate: bool = True
) -> bp.bibdatabase.BibDatabase:
    '''
    Returns a BibDatabase of all unique citations from a ProvDAG.
    If `deduplicate` is True references will be heuristically deduplicated.

    Parameters
    ----------
    dag : ProvDAG
        The ProvDAG object whose nodes contain citations to collect.
    deduplicate : bool
        Whether to deduplicate redundant citations.

    Returns
    -------
    bp.bibdatabase.BibDatabase
        A BibDatabase object containing the collected citations in bibtex
        format.
    '''
    bdb = bp.bibdatabase.BibDatabase()
    citations = []
    for node_uuid in dag:
        node = dag.get_node_data(node_uuid)
        # Skip no-prov nodes, which never have citations anyway
        if node is not None:
            node_citations = list(node.citations.values())
            citations.extend(node_citations)
    if deduplicate:
        citations = dedupe_citations(citations)

    bdb.entries = citations
    return bdb


class BibContent():
    '''
    A hashable data container capturing common bibtex fields

    Has many fields because keeping true duplicates is preferable to
    deduplicating two true non-duplicates.

    Parameters
    ----------
    entry : dict
        A dictionary of bibtex entries.
    '''
    def __init__(self, entry):
        self.title = entry.get('title'),
        self.author = entry.get('author'),
        self.journal = entry.get('journal')
        self.booktitle = entry.get('booktitle')
        self.year = entry.get('year')
        self.pages = entry.get('pages')

    def __eq__(self, other):
        return (
            type(self) == type(other) and
            self.title == other.title and
            self.author == other.author and
            self.journal == other.journal and
            self.booktitle == other.booktitle and
            self.year == other.year and
            self.pages == other.pages
        )

    def __hash__(self):
        return hash(
            str(self.title)
            + str(self.author)
            + str(self.journal)
            + str(self.journal)
            + str(self.booktitle)
            + str(self.year)
            + str(self.pages)
        )


def dedupe_citations(citations: List[Dict]) -> List[Dict]:
    '''
    Deduplicates citations based on bibtex id, bibtex content, and DOI.
    Citations are not guaranteed to be truly unique after deduplicating based
    on these values.

    Ensures only one qiime2 framework citation.

    Parameters
    ----------
    citations : list of dict
        The possibly redundant citations, each dict is a bibtex citation.

    Returns
    -------
    list of dict
        The deduplicated citations.
    '''
    deduped_citations = []
    is_framework_cited = False
    id_set = set()
    doi_set = set()
    content_set = set()
    for entry in citations:
        citation_id = entry['ID']

        if 'framework|qiime2' in citation_id:
            if not is_framework_cited:
                root = pkg_resources.resource_filename(
                    'qiime2.core.archive.provenance_lib', '.'
                )
                root = os.path.abspath(root)
                path = os.path.join(root, 'q2_citation.bib')
                with open(path) as bibtex_file:
                    q2_entry = bp.load(bibtex_file).entries.pop()

                q2_entry['ID'] = citation_id
                id_set.add(citation_id)
                deduped_citations.append(q2_entry)
                is_framework_cited = True
            continue

        # dedupe on id
        if citation_id in id_set:
            continue

        # dedupe on content
        entry_content = BibContent(entry)
        if entry_content in content_set:
            continue
        else:
            content_set.add(entry_content)

        # dedupe on doi if present
        doi = entry.get('doi')
        if doi is None:
            id_set.add(citation_id)
            deduped_citations.append(entry)
        elif doi not in doi_set:
            id_set.add(citation_id)
            doi_set.add(doi)
            deduped_citations.append(entry)

    return deduped_citations


def replay_citations(
    dag: ProvDAG,
    out_fp: str,
    deduplicate: bool = True,
    suppress_header: bool = False
):
    '''
    Writes a bibtex file containing all citations from a ProvDAG to disk.
    If `deduplicate` is True citations will be deduplicated, see
    `dedupe_citations()` for details.

    Parameters
    ----------
    dag : ProvDAG
        The provenance graph from which to collect citations.
    out_fp : str
        The filepath to which to write the citations.
    deduplicate : bool
        Whether to deduplicate the collected citations.
    suppress_header : bool
        Whether to forgo adding a header and footer to the output file.
    '''
    bib_db = collect_citations(dag, deduplicate=deduplicate)
    boundary = '#' * 79
    header = []
    footer = []
    extra = [
        '',
        '# This bibtex-formatted citation file can be imported into '
        'popular citation ',
        '# managers like Zotero and Mendeley, simplifying management and '
        'formatting.'
    ]
    if not suppress_header:
        header = build_header(boundary=boundary, extra_text=extra) + ['\n']
        footer = build_footer(dag=dag, boundary=boundary)
    if bib_db.entries_dict == {}:
        bib_db = 'No citations were registered to the used Actions.'
        with open(out_fp, 'w') as bibfile:
            bibfile.write(bib_db)
    else:
        with open(out_fp, 'w') as bibfile:
            bibfile.write('\n'.join(header))
            bibfile.write(BibTexWriter().write(bib_db))
            bibfile.write('\n'.join(footer))


def replay_supplement(
    payload: Union[str, ProvDAG],
    out_fp: str,
    validate_checksums: bool = True,
    parse_metadata: bool = True,
    use_recorded_metadata: bool = False,
    recurse: bool = False,
    deduplicate: bool = True,
    suppress_header: bool = False,
    verbose: bool = True,
    dump_recorded_metadata: bool = True
):
    '''
    Produces a zipfile package of useful documentation for in silico
    reproducibility of some QIIME 2 Result(s) from a ProvDAG, a QIIME 2
    Artifact, or a directory of Artifacts.

    Package includes:
    - replay scripts for all supported interfaces
    - a bibtex-formatted collection of all citations

    ProvDAG inputs retain their original config values. The
    `validate_checksums`, `parse_metadata`, `recurse`, and `verbose` parameters
    are disregarded if the payload is a ProvDAG.

    Parameters
    ----------
    payload : str or ProvDAG
        A filepath to an artifact or directory containing artifacts, or the
        ProvDAG to be parsed.
    out_fp : str
        The filepath to which to write the zip file.
    validate_checksums : bool
        Whether to perform checksum validation on the input artifact.
    parse_metadata : bool
        Whether to parse study metadata recorded in provenance.
    use_recorded_metadata : bool
        Whether to use the metadata recorded in provenance.
    recurse : bool
        Whether to recursively parse nested directories containing artifacts.
    deduplicate : bool
        Whether to deduplicate citations collected from provenance.
    suppress_header : bool
        Whether to forgo rendering the header and footer that are included
        by default in replay scripts.
    verbose : bool
        Whether to print status messages during processing.
    dump_recorded_metadata : bool
        Whether to write the metadata recorded in provenance to disk.
    '''
    dag = ProvDAG(
        artifact_data=payload,
        validate_checksums=validate_checksums,
        parse_metadata=parse_metadata,
        recurse=recurse,
        verbose=verbose
    )
    with tempfile.TemporaryDirectory() as tempdir:
        tempdir_path = pathlib.Path(tempdir)
        filenames = {
            'python3': 'python3_replay.py',
            'cli': 'cli_replay.sh',
        }

        for usage_driver in DRIVER_NAMES:
            md_out_fp = tempdir_path / 'recorded_metadata'
            rel_fp = filenames[usage_driver]
            tmp_fp = tempdir_path / rel_fp
            replay_provenance(
                payload=dag,
                out_fp=str(tmp_fp),
                usage_driver=usage_driver,
                use_recorded_metadata=use_recorded_metadata,
                suppress_header=suppress_header,
                verbose=verbose,
                dump_recorded_metadata=dump_recorded_metadata,
                md_out_fp=md_out_fp
            )
            print(f'The {usage_driver} replay script was written to {rel_fp}.')

        citations_fp = tempdir_path / 'citations.bib'
        replay_citations(
            dag,
            out_fp=str(citations_fp),
            deduplicate=deduplicate,
            suppress_header=suppress_header
        )
        print('The citations bibtex file was written to citations.bib.')

        out_fp = pathlib.Path(os.path.realpath(out_fp))
        if out_fp.suffix == '.zip':
            out_fp = out_fp.with_suffix('')

        shutil.make_archive(out_fp, 'zip', tempdir)
        print(f'The reproducibility package was written to {out_fp}.zip.')
