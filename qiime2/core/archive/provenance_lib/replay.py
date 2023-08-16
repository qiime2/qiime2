import bibtexparser as bp
from bibtexparser.bwriter import BibTexWriter
import networkx as nx
import os
import pathlib
import pkg_resources
import re
import shutil
import tempfile
from collections import UserDict
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Set, Union

from .archive_parser import ProvNode
from .parse import ProvDAG, UUID
from ._usage_drivers import (
    DRIVER_CHOICES, DRIVER_NAMES, SUPPORTED_USAGE_DRIVERS, Usage,
    build_header, build_footer
)
from .util import FileName, camel_to_snake
from ..provenance import MetadataInfo

from qiime2.sdk import PluginManager
from qiime2.sdk.usage import UsageVariable


@dataclass
class ReplayConfig():
    """
    fields:

    -use
      the usage driver to be used for provenance replay
    - dump_recorded_metadata
      If True, replay should write the metadata recorded in provenance to disk
      in .tsv format
    - use_recorded_metadata
      If True, replay should use the metadata recorded in provenance
    - pm
      an instance of the QIIME 2 PluginManager
    - md_context_has_been_printed
      a flag set by default and used internally, allows context to be printed
      once and only once.
    - no_provenance_context_has_been_printed
      indicates the no-provenance context documentation has not been printed
    - header
      if True, an introductory how-to header should be rendered to the script
    - verbose
      if True, progress will be reported to stdout
    - md_out_fp
      the directory path where caputred metadata should be written

    """
    use: Usage
    use_recorded_metadata: bool = False
    pm: PluginManager = PluginManager()
    md_context_has_been_printed: bool = False
    no_provenance_context_has_been_printed: bool = False
    header: bool = True
    verbose: bool = False
    dump_recorded_metadata: bool = True
    md_out_fp: FileName = ''


@dataclass
class ActionCollections():
    """
    std_actions are all normal, provenance-tracked q2 actions, arranged like:
    {<action_id>: {<output_node_uuid>: 'output_name',
                   <output_node_2_uuid:> 'output_name_2'},
     <action_2_id> : ...
     }

    no_provenance_nodes can't be organized by action, and in some cases we
    don't know anything but UUID for them, so we can fit what we need in a list
    """
    std_actions: Dict[UUID, Dict[UUID, str]] = field(default_factory=dict)
    no_provenance_nodes: List[UUID] = field(default_factory=list)


class UsageVarsDict(UserDict):
    """
    A dict where values are also unique. Used here as a UUID-queryable
    "namespace" of strings that can be passed to usage drivers for rendering
    into unique variable names.
    Non-unique values would cause namespace collisions.

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
    """
    def __setitem__(self, key: UUID, item: str) -> None:
        unique_item = self._uniquify(item)
        return super().__setitem__(key, unique_item)

    def _uniquify(self, var_name: str) -> str:
        """
        Appends _<some_int> to var_name, such that the returned name won't
        collide with any variable-name values that already exist in the dict.
        """
        some_int = 0
        unique_name = f"{var_name}_{some_int}"
        values = self.data.values()
        # no-prov nodes are stored with angle brackets around them, but
        # those brackets shouldn't be considered on uniqueness check
        while unique_name in values or f'<{unique_name}>' in values:
            some_int += 1
            unique_name = f"{var_name}_{some_int}"
        return unique_name

    def get_key(self, value: str):
        """
        Given some value in the dict, returns its key
        Results are predictable due to the uniqueness of dict values.

        Raises KeyError if search value does not exist.

        NOTE: If this proves too slow at scale, we can pivot to storing a
        second (reversed) dict for hashed lookups
        """
        for key, val in self.items():
            if value == val:
                return key
        raise KeyError(f"passed value '{value}' does not exist in this dict.")

    def wrap_val_in_angle_brackets(self, key: UUID):
        # TODO: Kinda unsafe. If we run it twice, it breaks uniquify.
        # Can we refactor this to behave like a repr instead of modifying the
        # value in the dict?
        super().__setitem__(key, f'<{self.data[key]}>')


@dataclass
class NamespaceCollections:
    usg_var_namespace: UsageVarsDict = field(default_factory=UsageVarsDict)
    usg_vars: Dict[UUID, UsageVariable] = field(default_factory=dict)
    action_namespace: Set[str] = field(default_factory=set)


def replay_provenance(payload: Union[FileName, ProvDAG],
                      out_fp: FileName,
                      usage_driver: DRIVER_CHOICES = 'python3',
                      validate_checksums: bool = True,
                      parse_metadata: bool = True,
                      recurse: bool = False,
                      use_recorded_metadata: bool = False,
                      suppress_header: bool = False,
                      verbose: bool = False,
                      dump_recorded_metadata: bool = True,
                      md_out_fp: FileName = '',):
    """
    Renders usage examples describing a ProvDAG, producing an interface-
    specific executable.

    Passed ProvDAGs retain their original config values.
    The following parameters are disregarded if payload is a ProvDAG,
    so may be left as default:
      - validate_checksums
      - parse_metadata
      - recurse
      - verbose
    """
    # Grab the right parse_metadata if the payload is already ProvDAG
    if hasattr(payload, 'cfg'):
        parse_metadata = payload.cfg.parse_study_metadata

    if use_recorded_metadata and not parse_metadata:
        raise ValueError(
            "Metadata not parsed for replay. Re-run with parse_metadata, or "
            "set use_recorded_metadata to False")

    if dump_recorded_metadata and not parse_metadata:
        raise ValueError(
            "Metadata not parsed, so cannot be written to disk. Re-run with "
            "parse_metadata, or set dump_recorded_metadata to False")

    if md_out_fp and not parse_metadata:
        raise ValueError(
            "Metadata not parsed, so cannot be written to disk. Re-run with "
            "parse_metadata, or do not pass a metadata output filepath "
            "argument.")

    if use_recorded_metadata and not dump_recorded_metadata:
        raise NotImplementedError(
            "In order to produce a replay script that uses metadata "
            "captured in provenance, that metadata must first be written to "
            "disk. Re-run with dump-recorded-metadata set to True, or "
            "use-recorded-metadata set to False. Possible future support for "
            "'touchless' replay from provenance is tracked in "
            "https://github.com/qiime2/provenance-lib/issues/98")

    # The ProvDAGParser handles ProvDAGs quickly, so we can just throw whatever
    # payload we get at this instead of maintaining per-data-type functions
    dag = ProvDAG(
        payload, validate_checksums, parse_metadata, recurse, verbose)

    cfg = ReplayConfig(use=SUPPORTED_USAGE_DRIVERS[usage_driver](),
                       use_recorded_metadata=use_recorded_metadata,
                       dump_recorded_metadata=dump_recorded_metadata,
                       verbose=verbose, md_out_fp=md_out_fp)
    # build order is handled by use.render, so doesn't matter here
    if not suppress_header:
        cfg.use.build_header()
        cfg.use.build_footer(dag)
    build_usage_examples(dag, cfg)

    if cfg.dump_recorded_metadata:
        print('metadata written to recorded_metadata/')

    output = cfg.use.render(flush=True)
    with open(out_fp, mode='w') as out_fh:
        out_fh.write(output)


def group_by_action(dag: ProvDAG, nodes: Iterator[UUID]) -> ActionCollections:
    """
    Provenance is organized around outputs, but replay cares about actions.
    This groups the nodes from a DAG by action, returning an ActionCollections
    aggregating the outputs related to each action.

    Takes an iterator of UUIDs, allowing us to influence the ordering of the
    grouping.

    In cases where a captured output_name is unavailable, we substitute the
    output data's Semantic Type, snake-cased because it will be used as
    a variable name if this data is rendered by ArtifactAPIUsage.
    """
    actions = ActionCollections()
    for node_id in nodes:
        if dag.node_has_provenance(node_id):
            data = dag.get_node_data(node_id)
            action_id = data.action._execution_details['uuid']

            output_name = data.action.output_name
            if output_name is None:
                output_name = camel_to_snake(data.type)

            try:
                actions.std_actions[action_id].update({node_id: output_name})
            except KeyError:
                actions.std_actions[action_id] = {node_id: output_name}
        else:
            actions.no_provenance_nodes.append(node_id)

    return actions


def build_usage_examples(dag: ProvDAG, cfg: ReplayConfig):
    """
    Builds a chained usage example representing the analysis `dag`.
    """
    usg_ns = NamespaceCollections()
    sorted_nodes = nx.topological_sort(dag.collapsed_view)
    actions = group_by_action(dag, sorted_nodes)

    for node_id in actions.no_provenance_nodes:
        n_data = dag.get_node_data(node_id)
        build_no_provenance_node_usage(n_data, node_id, usg_ns, cfg)

    for action_id in (std_actions := actions.std_actions):
        # We are replaying actions not nodes, so any associated node works
        some_node_id_from_this_action = next(iter(std_actions[action_id]))
        n_data = dag.get_node_data(some_node_id_from_this_action)
        if n_data.action.action_type == 'import':
            build_import_usage(n_data, usg_ns, cfg)
        else:
            build_action_usage(n_data, usg_ns, std_actions, action_id, cfg)


def build_no_provenance_node_usage(node: Optional[ProvNode],
                                   uuid: UUID,
                                   ns: NamespaceCollections,
                                   cfg: ReplayConfig):
    """
    Given a ProvNode (with no provenance), does something useful with it.
    Returns nothing, modifying the passed usage instance in place.

    # Basically:
    use.comment("Some context")
    use.comment("no-provenance nodes and descriptions")
    """
    if not cfg.no_provenance_context_has_been_printed:
        cfg.no_provenance_context_has_been_printed = True
        cfg.use.comment(
            "One or more nodes have no provenance, so full replay is "
            "impossible. Any commands we were able to reconstruct have been "
            "rendered, with the string descriptions below replacing actual "
            "inputs.")
        cfg.use.comment(
            "Original Node ID                       String Description")
    if node is None:
        # the node is a !no-provenance input and we have only UUID
        var_name = 'no-provenance-node'
    else:
        var_name = camel_to_snake(node.type)
    ns.usg_var_namespace.update({uuid: var_name})
    ns.usg_var_namespace.wrap_val_in_angle_brackets(uuid)

    # Make a usage variable for downstream consumption
    empty_var = cfg.use.usage_variable(
        ns.usg_var_namespace[uuid], lambda: None, 'artifact')
    ns.usg_vars.update({uuid: empty_var})

    # Log the no-prov node
    cfg.use.comment(f"{uuid}   {ns.usg_vars[uuid].to_interface_name()}")


def build_import_usage(node: ProvNode,
                       ns: NamespaceCollections,
                       cfg: ReplayConfig):
    """
    Given a ProvNode, adds an import usage example for it, roughly
    resembling the following.
    Returns nothing, modifying the passed usage instance in place.

    raw_seqs = use.init_format('raw_seqs', lambda: None, ext='fastq.gz')
    imported_seqs = use.import_from_format(
        'emp_single_end_sequences',
        'EMPSingleEndSequences',
        raw_seqs
    )

    The `lambda: None` is a placeholder for some actual data factory,
    and should not impact the rendered usage.
    """
    format_id = node._uuid + '_f'
    ns.usg_var_namespace.update({format_id: camel_to_snake(node.type) + '_f'})
    format_for_import = cfg.use.init_format(
        ns.usg_var_namespace[format_id], lambda: None)

    ns.usg_var_namespace.update({node._uuid: camel_to_snake(node.type)})
    use_var = cfg.use.import_from_format(
        ns.usg_var_namespace[node._uuid], node.type, format_for_import)
    ns.usg_vars.update({node._uuid: use_var})


def build_action_usage(node: ProvNode,
                       ns: NamespaceCollections,
                       std_actions: Dict[UUID, Dict[UUID, str]],
                       action_id: UUID,
                       cfg: ReplayConfig):
    """
    Adds an action usage example to `use` for some ProvNode.
    Returns nothing, modifying the passed usage instance in place.

    use.action(
        use.UsageAction(plugin_id='diversity_lib',
                        action_id='pielou_evenness'),
        use.UsageInputs(table=ft),
        use.UsageOutputNames(vector='pielou_vector')
    )
    """
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
                {unique_md_id: camel_to_snake(param_name)})
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
                    node, param_name, unique_md_id, ns.usg_var_namespace, cfg,
                    md_fn)
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
                    md = init_md_from_md_file(node, param_name, unique_md_id,
                                              ns.usg_var_namespace, cfg)
                else:
                    md = init_md_from_artifacts(param_val, ns, cfg)

            param_val = md

        inputs.update({param_name: param_val})

    usg_var = cfg.use.action(
        cfg.use.UsageAction(plugin_id=plugin, action_id=action),
        cfg.use.UsageInputs(**inputs),
        cfg.use.UsageOutputNames(**outputs))

    # write the usage vars into the UsageVars dict so we can use em downstream
    for res in usg_var:
        uuid_key = ns.usg_var_namespace.get_key(value=res.name)
        ns.usg_vars[uuid_key] = res


def _collect_action_inputs(ns: NamespaceCollections, node: ProvNode) -> dict:
    """
    Returns a dict containing the action Inputs from a ProvNode
    {input_name: input_vars}
    """
    inputs_dict = {}
    for input_name, uuids in node.action.inputs.items():
        # Some optional inputs take None as a default
        if uuids is not None:
            # Some inputs take collections of input strings, so:
            if type(uuids) is str:
                inputs_dict.update({input_name: ns.usg_vars[uuids]})
            else:  # it's a collection
                input_vars = []
                for uuid in uuids:
                    input_vars.append(ns.usg_vars[uuid])
                inputs_dict.update({input_name: input_vars})
    return inputs_dict


def _uniquify_output_names(ns: NamespaceCollections, raw_outputs) -> dict:
    """
    Returns a dict containing the uniquified output names from a ProvNode
    {output_name: uniquified_output_name}
    """
    outputs = {}
    for uuid, output_name in raw_outputs:
        ns.usg_var_namespace.update({uuid: output_name})
        uniquified_output_name = ns.usg_var_namespace[uuid]
        outputs.update({output_name: uniquified_output_name})
    return outputs


def init_md_from_recorded_md(node: ProvNode, param_name: str, md_id: str,
                             ns: UsageVarsDict, cfg: ReplayConfig,
                             md_fn: FileName) -> UsageVariable:
    """
    initializes and returns a Metadata UsageVariable with Metadata scraped
    and dumped to disk from provenance

    Assumes it will not be called if no metadata has been dumped to disk.
    This is enforced above in replay_provenance for faster failure

    Raises a ValueError if the node has no metadata
    """
    if not node.metadata:
        raise ValueError(
            'This function should only be called if the node has metadata.')
    parameter_name = ns[md_id][:-2]
    md_df = node.metadata[parameter_name]

    def factory():  # pragma: no cover
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


def init_md_from_md_file(node: ProvNode, param_name: str, md_id: str,
                         ns: UsageVarsDict, cfg: ReplayConfig) -> \
        UsageVariable:
    """
    initializes and returns a Metadata UsageVariable with no real data,
    mimicking a user passing md as a .tsv file
    """
    plugin = node.action.plugin
    action = node.action.action_name
    md = cfg.use.init_metadata(ns[md_id], lambda: None)
    if param_is_metadata_column(cfg, param_name, plugin, action):
        mdc_id = node._uuid + '_mdc'
        mdc_name = ns[md_id] + '_mdc'
        ns.update({mdc_id: mdc_name})
        md = cfg.use.get_metadata_column(ns[mdc_id], '<column name>', md)
    return md


def init_md_from_artifacts(md_inf: MetadataInfo,
                           ns: NamespaceCollections,
                           cfg: ReplayConfig) -> UsageVariable:
    """
    initializes and returns a Metadata UsageVariable with no real data,
    mimicking a user passing one or more QIIME 2 Artifacts as metadata

    We expect these usage vars are already in the namespace as artifacts if
    we're reading them in as metadata.
    """
    if not md_inf.input_artifact_uuids:
        raise ValueError("This funtion should not be used if"
                         "MetadataInfo.input_artifact_uuids is empty.")
    md_files_in = []
    for artif_id in md_inf.input_artifact_uuids:
        amd_id = artif_id + '_a'
        var_name = ns.usg_vars[artif_id].name + '_a'
        if amd_id not in ns.usg_var_namespace:
            ns.usg_var_namespace.update({amd_id: var_name})
            art_as_md = cfg.use.view_as_metadata(ns.usg_var_namespace[amd_id],
                                                 ns.usg_vars[artif_id])
            ns.usg_vars.update({amd_id: art_as_md})
        else:
            art_as_md = ns.usg_vars[amd_id]
        md_files_in.append(art_as_md)
    if len(md_inf.input_artifact_uuids) > 1:
        # We can't uniquify this normally, because one uuid can be merged with
        # combinations of others. One UUID does not a unique merge-id make.
        merge_id = '-'.join(md_inf.input_artifact_uuids)
        ns.usg_var_namespace.update({merge_id: 'merged_artifacts'})
        merged_md = cfg.use.merge_metadata(ns.usg_var_namespace[merge_id],
                                           *md_files_in)
        ns.usg_vars.update({merge_id: merged_md})
    return art_as_md


def dump_recorded_md_file(cfg: ReplayConfig, node: ProvNode, action_name: str,
                          md_id: str, fn: FileName):
    """
    Writes one metadata DataFrame passed to an action to .tsv
    Each action gets its own directory containing relevant md files.

    Raises a ValueError if the node has no metadata
    """
    # NOTE: node.metadata will also be None if no-parse-metadata is passed,
    # which would error unpredictably if the check preventing md dumping with
    # no-parse-provenance was removed from replay_provenance
    if node.metadata is None:
        raise ValueError(
            'This function should only be called if the node has metadata.')

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
        cfg: ReplayConfig, param: str, plg: str, action: str) -> bool:
    """
    Returns True if the param name `param` is registered as a MetadataColumn
    """
    try:
        plugin = cfg.pm.get_plugin(id=plg)
    except KeyError as e:
        msg = (re.sub("'", "", str(e)) +
               ' Visit library.qiime2.org to find plugins.')
        raise KeyError(msg)

    try:
        action_f = plugin.actions[action]
    except KeyError:
        raise KeyError(f'No action currently registered with id: {action}.')

    try:
        param_spec = action_f.signature.parameters[param]
    except KeyError:
        raise KeyError(f'No parameter currently registered with id: {param}')

    # HACK, but it works without relying on Q2's type system
    return ('MetadataColumn' in str(param_spec.qiime_type))


def uniquify_action_name(plugin: str, action: str, action_nmspace: set) -> str:
    """
    Creates a unique name by concatenating plugin_action_<counter>,
    and adds the name to action_namespace before returning it
    """
    counter = 0
    plg_action_name = f'{plugin}_{action}_{counter}'
    while plg_action_name in action_nmspace:
        counter += 1
        plg_action_name = f'{plugin}_{action}_{counter}'
    action_nmspace.add(plg_action_name)
    return plg_action_name


def collect_citations(dag: ProvDAG, deduplicate: bool = True) -> \
        bp.bibdatabase.BibDatabase:
    """
    Returns a BibDatabase of all unique citations from a ProvDAG.
    If `deduplicate`, refs will be heuristically deduplicated. e.g. by DOI
    """
    bdb = bp.bibdatabase.BibDatabase()
    cits = []
    for n_id in dag:
        p_node = dag.get_node_data(n_id)
        # Skip no-prov nodes, which never have citations anyway
        if p_node is not None:
            cit = list(p_node.citations.values())
            cits.extend(cit)
    if deduplicate:
        cits = dedupe_citations(cits)
    bdb.entries = cits
    return bdb


class BibContent():
    """
    A hashable data container capturing common bibtex fields

    Has many fields, b/c false negatives (i.e. self != other) are preferable.
    It is better to keep two dupes than to deduplicate a unique entry.
    """
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
            + str(self.pages))


def dedupe_citations(citations: List[Dict]) -> List[Dict]:
    """
    Heuristic attempts to reduce duplication in citations lists.
    E.g. capturing only one citation per DOI, ensuring one framework citation
    """
    dd_cits = []
    fw_cited = False
    id_set = set()
    doi_set = set()
    content_set = set()
    for entry in citations:
        id = entry['ID']
        # Deduplicate on bibtex key
        if id in id_set:
            continue

        # Write a single hardcoded framework citation
        if 'framework|qiime2' in id:
            if not fw_cited:
                root = pkg_resources.resource_filename(
                    'qiime2.core.archive.provenance_lib', '.')
                root = os.path.abspath(root)
                path = os.path.join(root, 'q2_citation.bib')
                with open(path) as bibtex_file:
                    q2_entry = bp.load(bibtex_file).entries.pop()

                q2_entry['ID'] = id
                id_set.add(id)
                dd_cits.append(q2_entry)
                fw_cited = True
            continue

        # deduplicate on content
        entry_content = BibContent(entry)
        if entry_content in content_set:
            continue
        else:
            content_set.add(entry_content)

        # Keep every unique entry without a doi
        if (doi := entry.get('doi')) is None:
            id_set.add(id)
            dd_cits.append(entry)
        # Keep one unique entry per non-framework doi
        else:
            if doi not in doi_set:
                id_set.add(id)
                dd_cits.append(entry)
                doi_set.add(doi)

    return dd_cits


def replay_citations(dag: ProvDAG, out_fp: FileName, deduplicate: bool = True,
                     suppress_header: bool = False):
    """
    Writes a .bib file representing all unique citations from a ProvDAG to disk
    If `deduplicate`, refs will be heuristically deduplicated. e.g. by DOI

    TODO: replay_citations_from_fp
    """
    bib_db = collect_citations(dag, deduplicate=deduplicate)
    boundary = '#' * 79
    header = []
    footer = []
    extra = [
        "",
        "# This bibtex-formatted citation file can be imported into "
        "popular citation ",
        "# managers like Zotero and Mendeley, simplifying management and "
        "formatting"
    ]
    if not suppress_header:
        header = build_header(boundary=boundary, extra_text=extra) + ['\n']
        footer = build_footer(dag=dag, boundary=boundary)
    if bib_db.entries_dict == {}:
        bib_db = "No citations were registered to the used Actions."
        with open(out_fp, 'w') as bibfile:
            bibfile.write(bib_db)
    else:
        with open(out_fp, 'w') as bibfile:
            bibfile.write('\n'.join(header))
            bibfile.write(BibTexWriter().write(bib_db))
            bibfile.write('\n'.join(footer))


def replay_supplement(payload: Union[FileName, ProvDAG],
                      out_fp: FileName,
                      validate_checksums: bool = True,
                      parse_metadata: bool = True,
                      use_recorded_metadata: bool = False,
                      recurse: bool = False,
                      deduplicate: bool = True,
                      suppress_header: bool = False,
                      verbose: bool = True,
                      dump_recorded_metadata: bool = True,
                      ):
    """
    Produces a zipfile package of useful documentation for in silico
    reproducibility of some QIIME 2 Result(s) from a ProvDAG, a QIIME 2
    Artifact, or a directory of Artifacts.

    Package includes:
    - replay scripts for all supported interfaces
    - a bibtex-formatted collection of all citations

    Passed ProvDAGs retain their original config values.
    The following parameters are disregarded if payload is a ProvDAG,
    so may be left as default:
      - validate_checksums
      - parse_metadata
      - recurse
    """
    # The ProvDAGParser handles ProvDAGs quickly, so we can just throw whatever
    # we get at this initializer instead of maintaining per-data-type functions
    dag = ProvDAG(artifact_data=payload, validate_checksums=validate_checksums,
                  parse_metadata=parse_metadata, recurse=recurse,
                  verbose=verbose)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = pathlib.Path(tmpdir)
        filenames = {
            'python3': 'python3_replay.py',
            'cli': 'cli_replay.sh',
        }

        for usage_driver in DRIVER_NAMES:
            md_out_fp = tmpdir_path / 'recorded_metadata'
            rel_fp = filenames[usage_driver]
            tmp_fp = tmpdir_path / rel_fp
            replay_provenance(
                payload=dag,
                out_fp=str(tmp_fp),
                usage_driver=usage_driver,
                use_recorded_metadata=use_recorded_metadata,
                suppress_header=suppress_header,
                verbose=verbose,
                dump_recorded_metadata=dump_recorded_metadata,
                md_out_fp=md_out_fp)
            print(f'{usage_driver} replay script written to {rel_fp}')

        tmp_fp = tmpdir_path / 'citations.bib'
        replay_citations(dag, out_fp=str(tmp_fp), deduplicate=deduplicate,
                         suppress_header=suppress_header)
        print('citations bibtex file written to citations.bib')

        out_fp = pathlib.Path(os.path.realpath(out_fp))
        # Drop .zip suffix if any so that we don't get some_file.zip.zip
        if out_fp.suffix == '.zip':
            out_fp = out_fp.with_suffix('')

        shutil.make_archive(out_fp, 'zip', tmpdir)
        print(f'reproducibility package written to {out_fp}.zip')
