# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from __future__ import annotations
import copy
from typing import Any, List, Optional, Set

import networkx as nx
import os
from pathlib import Path
from networkx.classes.reportviews import NodeView

from ._checksum_validator import ValidationCode, ChecksumDiff
from .archive_parser import (
    Config, ParserResults, ProvNode, Parser, ArchiveParser
)


class ProvDAG:
    '''
    A directed acyclic graph (DAG) representing the provenance of one or more
    QIIME 2 Archives (.qza or .qzv files).

    Parameters
    ----------
    artifact_data : Any
        The input data payload for the ProvDAG. Often the path to a file or
        directory on disk.
    validate_checksums : bool
        If True, Archives will be validated against their checksums.md5
        manifests.
    parse_metadata : bool
        If True, the metadata captured in the input Archives will be parsed
        and included in the ProvDAG.
    recurse : bool
        If True, and if artifact_data is a directory, will recursively parse
        all .qza and .qzv files within subdirectories.
    verbose : bool
        If True, will print parsed filenames to stdout, indicating progress.

    Attributes
    ----------
    dag : nx.DiGraph
        A Directed Acyclic Graph (DAG) representing the complete provenance
        of one or more QIIME 2 Artifacts. This DAG is comprehensive, including
        pipeline "alias" nodes as well as the inner nodes that compose each
        pipeline.
    parsed_artifact_uuids : Set[UUID]
        The set of user-passed terminal node uuids. Used to generate properties
        like `terminal_uuids`, this is a superset of terminal_uuids.
    terminal_uuids : Set[UUID]
        The set of terminal node ids present in the DAG, not including inner
        pipeline nodes.
    terminal_nodes : Set[ProvNode]
        The terminal ProvNodes present in the DAG, not including inner pipeline
        nodes.
    provenance_is_valid : ValidationCode
        The canonical indicator of provenance validity for a dag, this contains
        the lowest ValidationCode from all parsed Artifacts unioned into a
        given ProvDAG.
    checksum_diff : ChecksumDiff
        A ChecksumDiff representing all added, removed, and changed filepaths
        from all parsed Artifacts. If an artifact's checksums.md5 file is
        missing, this may be None. When multiple artifacts are unioned, this
        field prefers ChecksumDiffs over Nonetypes, which will be dropped.
    '''
    def __init__(
        self,
        artifact_data: Any = None,
        validate_checksums: bool = True,
        parse_metadata: bool = True,
        recurse: bool = False,
        verbose: bool = False,
    ):
        '''
        Creates a digraph by getting a parser from the parser dispatcher then
        parses the incoming data into a ParserResults, and then loads those
        results into key fields.
        '''
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
        return (
            'ProvDAG representing the provenance of the Artifacts: '
            f'{self._parsed_artifact_uuids}'
        )

    __str__ = __repr__

    def __len__(self) -> int:
        return len(self.dag)

    def __eq__(self, other) -> bool:
        if self.__class__ != other.__class__:
            return False
        if not nx.is_isomorphic(self.dag, other.dag):
            return False
        return True

    def __iter__(self):
        return iter(self.dag)

    @property
    def terminal_uuids(self) -> Set[str]:
        '''
        The UUIDs of the terminal nodes in the DAG, generated by selecting all
        nodes in a collapsed view of self.dag with an out-degree of zero.

        We memoize the set of terminal UUIDs to prevent unnecessary traversals,
        so must set self._terminal_uuid back to None in any method that
        modifies the structure of self.dag, or the nodes themselves (which are
        literal UUIDs). These methods include self.union() and
        self.relabel_nodes().
        '''
        if self._terminal_uuids is not None:
            return self._terminal_uuids
        cv = self.collapsed_view
        self._terminal_uuids = {
            uuid for uuid, out_degree in cv.out_degree()
            if out_degree == 0
        }
        return self._terminal_uuids

    @property
    def parsed_artifact_uuids(self) -> Set[str]:
        '''
        The set of user-passed terminal node uuids. Used to generate properties
        like self.terminal_uuids.
        '''
        return self._parsed_artifact_uuids

    @property
    def terminal_nodes(self) -> Set[ProvNode]:
        '''The terminal ProvNodes in the DAG's provenance.'''
        return {self.get_node_data(uuid) for uuid in self.terminal_uuids}

    @property
    def provenance_is_valid(self) -> ValidationCode:
        return self._provenance_is_valid

    @property
    def checksum_diff(self) -> Optional[ChecksumDiff]:
        return self._checksum_diff

    @property
    def nodes(self) -> NodeView:
        return self.dag.nodes

    @property
    def collapsed_view(self) -> nx.DiGraph:
        '''
        Returns a subsetted graphview of self.dag containing only nodes
        that are exterior to pipeline actions.
        '''
        outer_nodes = set()
        for terminal_uuid in self._parsed_artifact_uuids:
            outer_nodes |= self.get_outer_provenance_nodes(terminal_uuid)

        return nx.subgraph_view(self.dag, lambda node: node in outer_nodes)

    def has_edge(self, start_node: str, end_node: str) -> bool:
        return self.dag.has_edge(start_node, end_node)

    def node_has_provenance(self, uuid: str) -> bool:
        return self.dag.nodes[uuid]['has_provenance']

    def get_node_data(self, uuid: str) -> ProvNode:
        '''Returns a ProvNode from this ProvDAG selected by UUID.'''
        return self.dag.nodes[uuid]['node_data']

    def predecessors(self, node: str, dag: nx.DiGraph = None) -> Set[str]:
        '''
        Returns the parent UUIDs of a given node.

        Parameters
        ----------
        node : str
            The uuid of the node of interest.
        dag : nx.DiGraph
            The current provenance graph.

        Returns
        -------
        set of str
            The uuids of the parents of the node.
        '''
        dag = self.collapsed_view if dag is None else dag
        return set(self.dag.predecessors(node))

    @classmethod
    def union(cls, dags: List[ProvDAG]) -> ProvDAG:
        '''
        Class method that creates a new ProvDAG by unioning two or more graphs.

        The returned ProvDAG._parsed_artifact_uuids will include uuids from all
        dags, and other DAG attributes are reduced conservatively.

        Parameters
        ----------
        dags : list of ProvDAG
            A collection of ProvDAGs to union.

        Returns
        -------
        ProvDAG
            A single ProvDAG representing the union of all inputs.
        '''
        if len(dags) < 2:
            raise ValueError('Please pass at least two ProvDAGs.')

        union_dag = ProvDAG()
        union_dag.dag = nx.compose_all((dag.dag for dag in dags))

        union_dag._parsed_artifact_uuids = dags[0]._parsed_artifact_uuids
        union_dag._provenance_is_valid = dags[0]._provenance_is_valid
        union_dag._checksum_diff = dags[0].checksum_diff

        for next_dag in dags[1:]:
            union_dag._parsed_artifact_uuids = \
                union_dag._parsed_artifact_uuids.\
                union(next_dag._parsed_artifact_uuids)

            union_dag._provenance_is_valid = min(
                union_dag.provenance_is_valid, next_dag.provenance_is_valid
            )

            union_dag.cfg.parse_study_metadata = min(
                union_dag.cfg.parse_study_metadata,
                next_dag.cfg.parse_study_metadata
            )
            union_dag.cfg.perform_checksum_validation = min(
                union_dag.cfg.perform_checksum_validation,
                next_dag.cfg.perform_checksum_validation
            )

            # Here we retain as much data as possible, preferencing any
            # ChecksumDiff over None. This might mean we keep an empty
            # ChecksumDiff and drop a None ChecksumDiff, the latter of which is
            # used to indicate a missing checksums.md5 file in v5+ archives.
            # union_dag._provenance_is_valid will still be INVALID however.
            if next_dag.checksum_diff is None:
                continue

            if union_dag.checksum_diff is None:
                union_dag._checksum_diff = next_dag.checksum_diff
            else:
                union_dag.checksum_diff.added.update(
                    next_dag.checksum_diff.added
                )
                union_dag.checksum_diff.removed.update(
                    next_dag.checksum_diff.removed
                )
                union_dag.checksum_diff.changed.update(
                    next_dag.checksum_diff.changed
                )

        # make union._terminal_uuids be recalculated on next access
        union_dag._terminal_uuids = None
        return union_dag

    def get_outer_provenance_nodes(self, node_id: str = None) -> Set[str]:
        '''
        Performs depth-first traversal of a node's ancestors. Skips over nodes
        that are interior to a pipeline because pipeline output nodes point to
        the pipeline's inputs as parents not their direct parents inside
        the pipeline.

        Parameters
        ----------
        node_id : str
            The uuid of the node for which to discover ancestors.

        Returns
        -------
        set of str
            All ancestor uuids according to the above definition.
        '''
        nodes = set() if node_id is None else {node_id}
        parents = [edge_pair[0] for edge_pair in self.dag.in_edges(node_id)]
        for uuid in parents:
            nodes = nodes | self.get_outer_provenance_nodes(uuid)
        return nodes


# TODO: can this get nuked?
class EmptyParser(Parser):
    '''
    Creates empty ProvDAGs.
    Disregards Config, because it's not meaningful in this context.
    '''
    accepted_data_types = 'None'

    @classmethod
    def get_parser(cls, artifact_data: Any) -> Parser:
        if artifact_data is None:
            return EmptyParser()
        else:
            raise TypeError(f' in EmptyParser: {artifact_data} is not None.')

    def parse_prov(self, cfg: Config, data: None) -> ParserResults:
        '''
        Returns a static ParserResults with empty parsed_artifact_uuids,
        an empty graph, a valid ValidationCode, and a None ChecksumDiff.
        '''
        return ParserResults(
            parsed_artifact_uuids=set(),
            prov_digraph=nx.DiGraph(),
            provenance_is_valid=ValidationCode.VALID,
            checksum_diff=None,
        )


class DirectoryParser(Parser):
    accepted_data_types = \
        'filepath to a directory containing .qza/.qzv archives'

    @classmethod
    def get_parser(cls, artifact_data: Any) -> Parser:
        '''
        Return a DirectoryParser if appropriate.

        Parameters
        ----------
        artifact_data : Any
            Ideally a path to a directory containing one or more archives, but
            may be a different type during searches for other Parsers.

        Raises
        ------
        TypeError
            If something other than a str or path-like object is input.
        ValueError
            If the path does not point to a directory.
        '''
        try:
            is_dir = os.path.isdir(artifact_data)
        except TypeError:
            t = type(artifact_data)
            raise TypeError(
                f' in DirectoryParser: expects a directory, not a {t}.'
            )

        if not is_dir:
            raise ValueError(
                f' in DirectoryParser: {artifact_data} '
                'is not a valid directory.'
            )

        return DirectoryParser()

    def parse_prov(self, cfg: Config, data: str) -> ParserResults:
        '''
        Iterates over the directory's .qza and .qzv files, parsing them if
        their terminal node isn't already in the DAG.

        This behavior assumes that the ArchiveParsers capture all nodes within
        the archives they parse by default.

        Parameters
        ----------
        cfg : Config
            User-selected configuration options for whether to perform checksum
            validation, whether to parse study metadata, whether to recursively
            parse nested directories of artifacts, and whether to print status
            messages during processing.
        data : Any
            The path to the directory containing artifacts.

        Returns
        -------
        ParserResults
            A dataclass that stores the parsed artifact uuids, the parsed
            networkx graph, the provenance-is-valid flag, and the
            checksum diff.
        '''
        dir_name = Path(str(data).rstrip('/') + os.sep)
        if cfg.recurse:
            artifacts_to_parse = list(dir_name.rglob('*.qz[av]'))
            err_msg = (
                f'No .qza or .qzv files present in {dir_name} or any '
                'directory nested within it.'
            )
        else:
            artifacts_to_parse = list(dir_name.glob('*.qz[av]'))
            err_msg = (
                f'No .qza or .qzv files present in {dir_name}. Did you '
                'mean to recurse into nested directories?'
            )
        if not artifacts_to_parse:
            raise ValueError(err_msg)

        dag = ProvDAG()
        for archive in artifacts_to_parse:
            if cfg.verbose:
                print("parsing", archive)
            dag = ProvDAG.union([
                dag,
                ProvDAG(
                    archive,
                    cfg.perform_checksum_validation,
                    cfg.parse_study_metadata
                )
            ])

        return ParserResults(
            dag._parsed_artifact_uuids,
            dag.dag,
            dag.provenance_is_valid,
            dag.checksum_diff,
        )


def archive_not_parsed(root_uuid: str, dag: ProvDAG) -> bool:
    '''
    Checks if the archive with root_uuid has not already been parsed into dag
    which is defined as either not in the dag at all, or added only as a
    !no-provenance parent uuid.

    Parameters
    ----------
    root_uuid : str
        The root uuid of an archive.
    dag : ProvDAG
        The ProvDAG in which to search for the root uuid.

    Returns
    -------
    bool
        Indicating whether the archive represented by root_uuid has not been
        parsed into the ProvDAG.
    '''
    if root_uuid not in dag.dag:
        return True
    elif dag.get_node_data(root_uuid) is None:
        return True
    return False


class ProvDAGParser(Parser):
    '''
    Effectively a ProvDAG copy constructor, this "parses" a ProvDAG, loading
    its data into a new ProvDAG.

    Disregards Config, because it's not meaningful in this context.
    '''
    accepted_data_types = 'ProvDAG'

    @classmethod
    def get_parser(cls, artifact_data: Any) -> Parser:
        '''
        Returns ProvDAGParser if appropriate.

        Parameters
        ----------
        artifact_data : Any
            Hopefully a ProvDAG but may be a different type during searches for
            the proper Parser.

        Returns
        -------
        ProvDAGParser
            An instance of ProvDAGParser if artifact_data is a ProvDAG.

        Raises
        ------
        TypeError
            If artifact_data is not a ProvDAG.
        '''
        if isinstance(artifact_data, ProvDAG):
            return ProvDAGParser()
        else:
            raise TypeError(
                f' in ProvDAGParser: {artifact_data} is not a ProvDAG.'
            )

    def parse_prov(self, cfg: Config, dag: ProvDAG) -> ParserResults:
        '''
        Parses a ProvDAG returning a ParserResults by deep copying existing
        attributes that live on the ProvDAG and make up a ParserResults.

        Parameters
        ----------
        cfg : Config
            Ignored because a ProvDAG is not being constructed from scratch.
            Present for inheritance purposes.
        dag : ProvDAG
            The ProvDAG to parse, read: copy attributes from.

        Returns
        -------
        ParserResults
            A dataclass that stores the parsed artifact uuids, the parsed
            networkx graph, the provenance-is-valid flag, and the
            checksum diff.
        '''
        return ParserResults(
            copy.deepcopy(dag._parsed_artifact_uuids),
            copy.deepcopy(dag.dag),
            copy.deepcopy(dag.provenance_is_valid),
            copy.deepcopy(dag.checksum_diff),
        )


def parse_provenance(cfg: Config, payload: Any) -> ParserResults:
    '''
    Parses some data payload into a ParserResults object ingestible by ProvDAG.

    Parameters
    ----------
    cfg : Config
        A dataclass that stores four boolean flags: whether to perform
        checksum validation, whether to parse study metadata, whether to
        recursively parse nested directories, and whether to enable verbose
        mode.
    payload : Any
        The payload to attempt to parse, commonly a path to an archive or
        directory containing archives.

    Returns
    -------
    ParserResults
        A dataclass that stores the parsed artifact uuids, the parsed
        networkx graph, the provenance-is-valid flag, and the
        checksum diff.

    '''
    parser = select_parser(payload)
    return parser.parse_prov(cfg, payload)


def select_parser(payload: Any) -> Parser:
    '''
    Attempts to find a parser that can handle some given payload.

    Parameters
    ----------
    payload : Any
        The payload for which to find a parser.

    Returns
    -------
    Parser
        The appropriate Parser for the payload type.

    Raises
    ------
    UnparseableDataError
        If no appropriate parser could be found for the payload.
    '''
    _PARSER_TYPE_REGISTRY = [
        ArchiveParser,
        DirectoryParser,
        ProvDAGParser,
        EmptyParser
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
