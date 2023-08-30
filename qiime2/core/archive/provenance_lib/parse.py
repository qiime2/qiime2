from __future__ import annotations
import copy
from typing import Any, List, Mapping, Optional, Set

# import glob
import networkx as nx
import os
from pathlib import Path
from networkx.classes.reportviews import NodeView
import zipfile

from . import _checksum_validator
from .archive_parser import (
    Config, ParserResults, ProvNode, Parser, ArchiveParser
)
from .util import get_root_uuid


class ProvDAG:
    """
    A directed acyclic graph (DAG) representing the provenance of one or more
    QIIME 2 Archives (.qza or .qzv files). It is safest to manipulate ProvDAGs
    with the properties and methods exposed here, but advanced users may
    get value from using native methods of the core networkx.DiGraph at `.dag`.

    ## Parameters

    artifact_data: Any = None - the input data payload for the ProvDAG. Often
        the path to a file or directory on disk
    validate_checksums: bool = True - if True, Archives will be validated
        against their checksums.md5 manifests
    parse_metadata: bool = True - if True, the metadata captured in the input
        Archives will be parsed and included in the ProvDAG
    recurse: bool = False - if True, and if artifact_data is a directory,
        will recursively parse all .qza and .qzv files within subdirectories
    verbose: bool = False - if True, will print parsed filenames to stdout,
        indicating progress

    ## Properties

    dag: nx.DiGraph - a Directed Acyclic Graph (DAG) representing the complete
        provenance of one or more QIIME 2 Artifacts. This DAG is comprehensive,
        including pipeline "alias" nodes as well as the inner nodes that
        compose each pipeline.
    parsed_artifact_uuids: Set[UUID] - the set of user-passed terminal node
        uuids. Used to generate properties like `terminal_uuids`, this is a
        superset of terminal_uuids.
    terminal_uuids: Set[UUID] - the set of terminal node ids present in the
        DAG, not including inner pipeline nodes.
    terminal_nodes: Set[ProvNode] - the terminal ProvNodes present in the DAG,
        not including inner pipeline nodes.
    provenance_is_valid: checksum_validator.ValidationCode - the canonical
        indicator of provenance validity for a dag, this contains the _poorest_
        ValidationCode from all parsed Artifacts unioned into a given ProvDAG.
    checksum_diff: checksum_validator.ChecksumDiff - a ChecksumDiff
        representing all added, removed, and changed filepaths from all parsed
        Artifacts. If an artifact's checksums.md5 file is missing, this may
        be None. When multiple artifacts are unioned, this field prefers
        ChecksumDiffs over Nonetypes, which will be dropped. For this reason,
        provenance_is_valid is a more reliable indicator of provenance validity
        than checksum_diff.

    ## Methods

    nodes: networkx.classes.reportview.NodeView - A NodeView of self.dag
    relabel_nodes: Optional[ProvDAG]: takes an {old_node: new_node} mapping,
        and relabels the nodes in self.dag.
    union: ProvDAG - a class method that returns the union of multiple ProvDAGs

    ## Builtin support

    len: int - ProvDAG supports the builtin `len` just as nx.DiGraph does,
        with `len(mydag)` returning the number of nodes in `mydag.dag`
    __eq__: bool - checks that self and other are of the same class and
        have isomorphic DAGs
    __iter__: iterator over the underlying nx.DiGraph (ProvDAG.dag)

    ## Available GraphViews
    Graphviews are read-only subgraphs of networkx graphs. They behave just
    like DiGraphs but are cheaper. If you make many views of views, they lag.

    complete: `mydag.dag` is the DiGraph containing all recorded provenance
        nodes for this ProvDAG
    collapsed_view: `mydag.collapsed_view` returns a DiGraph (GraphView)
    containing a node for each standalone Action or Visualizer and one single
    node for each Pipeline (like q2view provenance trees)

    ## About the Nodes

    DiGraph nodes are literally UUIDs (strings)

    Every node has the following attributes:
    node_data: Optional[ProvNode]
    has_provenance: bool

    ## No-provenance nodes

    NOTE: The following is only relevant if some of your provenance is from a
    QIIME 2 version earlier than v2.0.6, which was superseded on 10/24/16.

    When parsing v1+ archives, v0 ancestor nodes without tracked provenance
    (e.g. !no-provenance inputs) are discovered only as parents to the current
    inputs. They are added to the DAG when we add in-edges to "real" provenance
    nodes. These nodes are explicitly assigned the node attributes above,
    allowing red-flagging of no-provenance nodes, as all nodes have a
    has_provenance attribute.

    A node with has_provenance = False may be a "complete" v0 node or a
    uuid-only !no-provenance input.

    A node with node_data = None will always be a uuid-only !no-provence input

    No-provenance nodes with no v1+ children will always appear as disconnected
    members of the DiGraph.
    """
    def __init__(
        self,
        artifact_data: Any = None,
        validate_checksums: bool = True,
        parse_metadata: bool = True,
        recurse: bool = False,
        verbose: bool = False,
    ):
        """
        Creates a ProvDAG (digraph) by getting a parser from the parser
        dispatcher, using it to parse the incoming data into a ParserResults,
        and then loading those Results into key fields.
        """
        cfg = Config(validate_checksums, parse_metadata, recurse, verbose)
        parser_results = parse_provenance(cfg, artifact_data)

        self.cfg = cfg
        self._parsed_artifact_uuids = parser_results.parsed_artifact_uuids
        self.dag = parser_results.prov_digraph
        self._provenance_is_valid = parser_results.provenance_is_valid
        self._checksum_diff = parser_results.checksum_diff

        # clear cache whenever we create a new ProvDAG
        self._terminal_uuids = None

    def __repr__(self) -> str:
        return ('ProvDAG representing these Artifacts '
                f'{self._parsed_artifact_uuids}')

    __str__ = __repr__

    def __len__(self) -> int:
        return len(self.dag)

    def __eq__(self, other) -> bool:
        if (self.__class__ != other.__class__ or
                not nx.is_isomorphic(self.dag, other.dag)):
            return False
        else:
            return True

    def __iter__(self):
        return iter(self.dag)

    @property
    def terminal_uuids(self) -> Set[str]:
        """
        The UUIDs of the terminal nodes in the DAG, generated by selecting all
        nodes in a collapsed view of self.dag with an out-degree of zero.

        We memoize the set of terminal UUIDs to prevent unnecessary traversals,
        so must set self._terminal_uuid back to None in any method that
        modifies the structure of self.dag, or the nodes themselves (which are
        literal UUIDs). These methods include at least union and relabel_nodes.
        """
        if self._terminal_uuids is not None:
            return self._terminal_uuids
        cv = self.collapsed_view
        self._terminal_uuids = {uuid for uuid, out_degree in cv.out_degree()
                                if out_degree == 0}
        return self._terminal_uuids

    @property
    def parsed_artifact_uuids(self) -> Set[str]:
        """
        The the set of user-passed terminal node uuids. Used to generate
        properties like `terminal_uuids`, this is a superset of
        terminal_uuids."""
        return self._parsed_artifact_uuids

    @property
    def terminal_nodes(self) -> Set[ProvNode]:
        """The terminal ProvNodes in the DAG's provenance"""
        return {self.get_node_data(uuid) for uuid in self.terminal_uuids}

    @property
    def provenance_is_valid(self) -> _checksum_validator.ValidationCode:
        return self._provenance_is_valid

    @property
    def checksum_diff(self) -> Optional[_checksum_validator.ChecksumDiff]:
        return self._checksum_diff

    @property
    def nodes(self) -> NodeView:
        return self.dag.nodes

    @property
    # NOTE: This actually returns a graphview, which is a read-only DiGraph
    def collapsed_view(self) -> nx.DiGraph:
        outer_nodes = set()
        for terminal_uuid in self._parsed_artifact_uuids:
            outer_nodes |= self.get_outer_provenance_nodes(terminal_uuid)

        def n_filter(node):
            return node in outer_nodes

        return nx.subgraph_view(self.dag, filter_node=n_filter)

    def has_edge(self, start_node: str, end_node: str) -> bool:
        """
        Returns True if the edge u, v is in the graph
        Calls nx.DiGraph.has_edge
        """
        return self.dag.has_edge(start_node, end_node)

    def node_has_provenance(self, uuid: str) -> bool:
        return self.dag.nodes[uuid]['has_provenance']

    def get_node_data(self, uuid: str) -> ProvNode:
        """ Returns a ProvNode from this ProvDAG selected by UUID """
        return self.dag.nodes[uuid]['node_data']

    def predecessors(self, node: str, dag: nx.DiGraph = None) \
            -> Set[str]:
        """ Returns the parent UUIDs of a given node """
        dag = self.collapsed_view if dag is None else dag
        return set(self.dag.predecessors(node))

    def relabel_nodes(self, mapping: Mapping, copy: bool = False) -> \
            Optional[ProvDAG]:
        """
        Helper method for safe use of nx.relabel.relabel_nodes.
        By default, this updates the labels of self.dag in place.
        With copy=True, returns a copy of self with nodes relabeled.

        Also updates the DAG's _parsed_artifact_uuids to match the new labels
        to head off KeyErrors downstream, and clears the _terminal_uuids cache.
        """
        mod_dag = self
        if copy:
            mod_dag = ProvDAG(self)

        # rename node uuids in the provnode data payloads for consistency
        for node_id in mod_dag:
            mod_dag.get_node_data(node_id)._uuid = mapping[node_id]

        # then update the dag itself
        nx.relabel_nodes(mod_dag.dag, mapping, copy=False)

        mod_dag._parsed_artifact_uuids = {mapping[uuid] for
                                          uuid in self._parsed_artifact_uuids}

        # Clear the _terminal_uuids cache of the dag whose nodes we're changing
        # so that property returns correctly
        mod_dag._terminal_uuids = None

        if copy:
            return mod_dag
        else:
            return None

    @classmethod
    def union(cls, dags: List[ProvDAG]) -> ProvDAG:
        """
        Class method that creates a new ProvDAG by unioning the graphs in an
        arbitrary number of ProvDAGs.

        The returned DAG's _parsed_artifact_uuids will include uuids from all
        dags, and other DAG attributes are reduced conservatively.
        """
        if len(dags) < 2:
            raise ValueError("Please pass at least two ProvDAGs")

        new_dag = ProvDAG()
        new_dag.dag = nx.compose_all((dag.dag for dag in dags))

        new_dag._parsed_artifact_uuids = dags[0]._parsed_artifact_uuids
        new_dag._provenance_is_valid = dags[0]._provenance_is_valid
        new_dag._checksum_diff = dags[0].checksum_diff

        for dag in dags[1:]:
            new_dag._parsed_artifact_uuids = new_dag._parsed_artifact_uuids \
                                             .union(dag._parsed_artifact_uuids)

            # provenance of union is only valid if provenances of all members
            # are valid
            new_dag._provenance_is_valid = min(new_dag.provenance_is_valid,
                                               dag.provenance_is_valid)

            new_dag.cfg.parse_study_metadata = min(
                new_dag.cfg.parse_study_metadata,
                dag.cfg.parse_study_metadata)
            new_dag.cfg.perform_checksum_validation = min(
                new_dag.cfg.perform_checksum_validation,
                dag.cfg.perform_checksum_validation)

            # Here we retain as much data as possible, preferencing any
            # ChecksumDiff over None. This might mean we keep a clean/empty
            # ChecksumDiff and drop None, used to indicate a missing
            # checksums.md5 file in a v5+ archive. _provenance_is_valid will
            # still be INVALID in this case.
            if dag.checksum_diff is None:
                continue

            if new_dag.checksum_diff is None:
                new_dag._checksum_diff = dag.checksum_diff
            else:
                new_dag.checksum_diff.added.update(dag.checksum_diff.added)
                new_dag.checksum_diff.removed.update(dag.checksum_diff.removed)
                new_dag.checksum_diff.changed.update(dag.checksum_diff.changed)

        # Clear the _terminal_uuids cache so that property returns correctly
        new_dag._terminal_uuids = None
        return new_dag

    def get_outer_provenance_nodes(self, _node_id: str = None) -> Set[str]:
        """
        Selective depth-first traversal of this node_id's ancestors.
        Returns the set of "outer" nodes that represent "collapsed" provenance
        like that seen in q2view (i.e. all standalone Actions and Visualizers,
        and a single node for each Pipeline).

        Because the terminal/alias nodes created by pipelines show _pipeline_
        inputs, this recursion skips over all inner nodes.

        NOTE: _node_id exists to support recursive calls and may produce
        unexpected results if e.g. an "inner" node ID is passed by an external
        caller.
        """
        nodes = set() if _node_id is None else {_node_id}
        parents = [edge_pair[0] for edge_pair in self.dag.in_edges(_node_id)]
        for uuid in parents:
            nodes = nodes | self.get_outer_provenance_nodes(uuid)
        return nodes


class EmptyParser(Parser):
    """
    Creates empty ProvDAGs.
    Disregards Config, because it's not meaningful in this context.
    """
    accepted_data_types = "None"

    @classmethod
    def get_parser(cls, artifact_data: Any) -> Parser:
        if artifact_data is None:
            return EmptyParser()
        else:
            raise TypeError(f" in EmptyParser: {artifact_data} is not None")

    def parse_prov(self, cfg: Config, data: None) -> ParserResults:
        return ParserResults(
            parsed_artifact_uuids=set(),
            prov_digraph=nx.DiGraph(),
            provenance_is_valid=_checksum_validator.ValidationCode.VALID,
            checksum_diff=None,
        )


class DirectoryParser(Parser):
    accepted_data_types = \
        'filepath to a directory containing .qza/.qzv archives'

    @classmethod
    def get_parser(cls, artifact_data: Any) -> Parser:
        """
        Return the appropriate Parser if this Parser type can handle the data
        passed in.

        Should raise an appropriate exception if this Parser cannot handle the
        data.
        """
        try:
            is_dir = os.path.isdir(artifact_data)
        except TypeError:
            t = type(artifact_data)
            raise ValueError(
                f" in DirectoryParser: expects a directory, not a {t}")

        if not is_dir:
            raise ValueError(f" in DirectoryParser: {artifact_data} "
                             "is not a valid directory.")

        return DirectoryParser()

    def parse_prov(self, cfg: Config, data: Any) -> ParserResults:
        """
        Iterates over the directory's .qza and .qzv files, parsing them if
        their terminal node isn't already in the DAG.

        This behavior assumes that the ArchiveParsers capture all nodes
        within the archives they parse by default.
        """
        dir_name = Path(str(data).rstrip('/') + os.sep)
        if cfg.recurse:
            # "empty" generators don't fail if not checks, so cast to list
            artifacts_to_parse = list(dir_name.rglob('*.qz[av]'))
        else:
            artifacts_to_parse = list(dir_name.glob('*.qz[av]'))
        if not artifacts_to_parse:
            raise ValueError(f"No .qza or .qzv files present in {dir_name}")

        dag = ProvDAG()
        for archive in artifacts_to_parse:
            if cfg.verbose:
                print("parsing", archive)
            with zipfile.ZipFile(archive) as zf:
                root_id = get_root_uuid(zf)
            if archive_not_parsed(root_id, dag):
                dag = ProvDAG.union(
                    [dag,
                     ProvDAG(archive,
                             cfg.perform_checksum_validation,
                             cfg.parse_study_metadata)])
            else:
                # Even if we skip a redundant file for efficiency,
                # we should add its UUID to the list of parsed artifacts
                dag._parsed_artifact_uuids.add(root_id)

        return ParserResults(
            dag._parsed_artifact_uuids,
            dag.dag,
            dag.provenance_is_valid,
            dag.checksum_diff,
            )


def archive_not_parsed(root_id: str, dag: ProvDAG) -> bool:
    """
    returns True if the archive with root_uuid has not been parsed into dag
    (not in the dag at all, or added only as a !no-provenance parent id)
    """
    return (root_id in dag.dag and dag.get_node_data(root_id) is None) \
        or root_id not in dag.dag


class ProvDAGParser(Parser):
    """
    Effectively a ProvDAG copy constructor, this "parses" a ProvDAG, loading
    its data into a new ProvDAG.

    Disregards Config, because it's not meaningful in this context.
    """
    accepted_data_types = 'ProvDAG'

    @classmethod
    def get_parser(cls, artifact_data: Any) -> Parser:
        if isinstance(artifact_data, ProvDAG):
            return ProvDAGParser()
        else:
            raise TypeError(
                f" in ProvDAGParser: {artifact_data} is not a ProvDAG")

    def parse_prov(self, cfg: Config, pdag: ProvDAG) -> ParserResults:
        return ParserResults(
            copy.deepcopy(pdag._parsed_artifact_uuids),
            copy.deepcopy(pdag.dag),
            copy.deepcopy(pdag.provenance_is_valid),
            copy.deepcopy(pdag.checksum_diff),
        )


def parse_provenance(cfg: Config, payload: Any) -> ParserResults:
    """
    Parses some data payload into a ParserResults object ingestible by ProvDAG.
    """
    parser = select_parser(payload)
    return parser.parse_prov(cfg, payload)


def select_parser(payload: Any) -> Parser:
    '''
    Attempts to find a parser that can handle some given payload. Raises
    '''
    _PARSER_TYPE_REGISTRY = [
        ArchiveParser,
        DirectoryParser,
        ProvDAGParser,
        EmptyParser,
    ]

    accepted_data_types = [
        parser.accepted_data_types for parser in _PARSER_TYPE_REGISTRY
    ]

    optional_parser = None
    errors = []
    for parser in _PARSER_TYPE_REGISTRY:
        try:
            optional_parser = parser.get_parser(payload)
            if optional_parser is not None:
                return optional_parser
        except Exception as e:
            errors.append(e)

    err_msg = (
        f'Input data {payload} is not supported.\n'
        'Parsers are available for the following data types: '
        f'{accepted_data_types}.\n'
        'The following errors were caught while trying to identify a parser '
        'that can_handle this input data:\n'
    )
    for e in errors:
        err_msg += str(type(e)) + str(e) + '\n'
    raise UnparseableDataError(err_msg)


class UnparseableDataError(Exception):
    pass