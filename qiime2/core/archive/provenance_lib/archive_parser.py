import abc
import os
import pandas as pd
import pathlib
import tempfile
import yaml
import warnings
from zipfile import ZipFile

from dataclasses import dataclass
from datetime import timedelta
from io import BytesIO
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import bibtexparser as bp
import networkx as nx

from ._checksum_validator import (
    ValidationCode, ChecksumDiff, validate_checksums
)
from .util import get_root_uuid, get_nonroot_uuid, parse_version
from ..provenance import MetadataInfo


@dataclass
class Config():
    '''
    Dataclass that stores user-selected configuration options.

    Attributes
    ----------
    perform_checksum_validation : bool
        Whether to opt in or out of checksum validation.
    parse_study_metadata : bool
        Whether to parse study metadata stored in provenance.
    recurse : bool
        Whether to recursively parse nested directories that contain artifacts.
    verbose : bool
        Whether to print status messages to stdout during processing.
    '''
    perform_checksum_validation: bool = True
    parse_study_metadata: bool = True
    recurse: bool = False
    verbose: bool = False


@dataclass
class ParserResults():
    '''
    Results generated and returned by a ParserVx.

    Attributes
    ----------
    parsed_artifact_uuids : set of str
        The uuids of the artifacts directly parsed by a parser. Does not
        include the uuids of artifact parsed from provenance. When parsing
        a single archive this is a single member set of that uuid. When
        parsing a directory, it is the set of all artifact uuids in that
        directory.
    prov_digraph : nx.Digraph
        The directed acyclic graph representation of the parsed provenance as
        an nx.DiGraph object.
    provenance_is_valid : ValidationCode
        A flag indicating the level of checksum validation.
    checksum_diff : ChecksumDiff or None
        A tuple of three dictionaries indicating the uuids of files that have
        been 1) added 2) removed or 3) changed in the archive since the
        archive was checksummed.
        None if no checksum validation was perfomed, e.g. when opted out or
        impossible because archive version did not support checksums, or when
        checksums.md5 missing from archive where it was expected.
        Interpretable only in conjunction with provenance_is_valid.
    '''
    parsed_artifact_uuids: Set[str]
    prov_digraph: nx.DiGraph
    provenance_is_valid: ValidationCode
    checksum_diff: Optional[ChecksumDiff]


class ProvNode:
    '''
    One node of a provenance DAG, describing one QIIME2 Result.
    '''

    @property
    def _uuid(self) -> str:
        return self._result_md.uuid

    @_uuid.setter
    def _uuid(self, new_uuid: str):
        '''
        ProvNode's UUID. Safe for use as getter. Prefer ProvDAG.relabel_nodes
        as a setter because it preserves alignment between ids across the dag
        and its ProvNodes.
        '''
        self._result_md.uuid = new_uuid

    @property
    def type(self) -> str:
        return self._result_md.type

    @property
    def format(self) -> Optional[str]:
        return self._result_md.format

    @property
    def archive_version(self) -> str:
        return self._archive_version

    @property
    def framework_version(self) -> str:
        return self._framework_version

    @property
    def has_provenance(self) -> bool:
        return int(self.archive_version) > 1

    @property
    def citations(self) -> Dict:
        citations = {}
        if hasattr(self, '_citations'):
            citations = self._citations.citations
        return citations

    @property
    def metadata(self) -> Optional[Dict[str, pd.DataFrame]]:
        '''
        A dict containing {parameter_name: metadata_dataframe} pairs where
        parameter_name is the registered name of the parameter the Metadata
        or MetadataColumn was passed to.

        Returns an empty dict if this action takes no Metadata or
        MetadataColumn.

        Returns None if this action has no metadata because the archive has no
        provenance, or the user opted out of metadata parsing.
        '''
        self._metadata: Optional[Dict[str, pd.DataFrame]]

        md = None
        if hasattr(self, '_metadata'):
            md = self._metadata
        return md

    @property
    def _parents(self) -> Optional[List[Dict[str, str]]]:
        '''
        A list of single-item {Type: UUID} dicts describing this
        action's inputs, including Artifacts passed as Metadata parameters.

        Returns [] if this action is an Import.

        NOTE: This property is private because it is slightly unsafe,
        reporting original node IDs that are not updated if the user renames
        nodes using the networkx API instead of ProvDAG.relabel_nodes.
        ProvDAG and its extensions should use the networkx.DiGraph itself to
        work with ancestry when possible.
        '''
        if not self.has_provenance:
            return None

        inputs = self.action._action_details.get('inputs')
        parents = []
        if inputs is not None:
            # Inputs are a list of single-item dicts
            for input in inputs:
                (name, value), = input.items()
                # value is usually a uuid, but may be a collection of uuids
                # the following are specced in qiime2/core/type/collection
                if type(value) in (set, list, tuple):
                    for i in range(len(value)):
                        # Make these unique in case the single-item dicts get
                        # merged into a single dict downstream
                        if type(value[i]) is dict:
                            unq_name, = value[i].keys()
                            v, = value[i].values()
                        else:
                            unq_name = f'{name}_{i}'
                            v = value[i]
                        parents.append({unq_name: v})
                elif value is not None:
                    parents.append({name: value})
                else:
                    # skip None-by-default optional inputs
                    pass

        return parents + self._artifacts_passed_as_md

    def __init__(
        self,
        cfg: Config,
        zf: ZipFile,
        node_fps: List[pathlib.Path]
    ):
        '''
        Constructs a ProvNode from a zipfile and the collected
        provenance-relevant filepaths for a single result within it.
        '''
        for fp in node_fps:
            if fp.name == 'VERSION':
                self._archive_version, self._framework_version = \
                    parse_version(zf)
            elif fp.name == 'metadata.yaml':
                self._result_md = _ResultMetadata(zf, str(fp))
            elif fp.name == 'action.yaml':
                self.action = _Action(zf, str(fp))
            elif fp.name == 'citations.bib':
                self._citations = _Citations(zf, str(fp))
            elif fp.name == 'checksums.md5':
                # Handled in ProvDAG
                pass

        if self.has_provenance:
            all_metadata_fps, self._artifacts_passed_as_md = \
                self._get_metadata_from_Action(self.action._action_details)
            if cfg.parse_study_metadata:
                self._metadata = self._parse_metadata(zf, all_metadata_fps)

    def _get_metadata_from_Action(
        self, action_details: Dict[str, List]
    ) -> Tuple[Dict[str, str], List[Dict[str, str]]]:
        '''
        Gathers data related to Metadata and MetadataColumn-based metadata
        files from the parsed action.yaml file.

        Captures filepath and parameter-name data for all study metadata
        files, so that these can be located for parsing, and then associated
        with the correct parameters during replay. It captures uuids for all
        artifacts passed to this action as metadata so they can be included as
        parents of this node.

        Parameters
        ----------
        action_details : dict
            The parsed dictionary of the `action` section from action.yaml.

        Returns
        -------
        tuple of (all_metadata, artifacts_as_metadata)
            Where all_metadata is a dict of
            {parameter_name: filename}.
            Where artifacts_as_metadata is a list of single-items dict of the
            structure {'artifact_passed_as_metadata': <uuid>}.

        Notes
        -----
        When Artifacts are passed as Metadata, they are captured in
        action['parameters'], rather than in action['inputs'] with the other
        Artifacts. Semantic Type data is thus not captured. This function
        returns a filler 'Type' for all UUIDs discovered here:
        'artifact_passed_as_metadata'. Because Artifacts passed (viewed) as
        Metadata retain their provenance, downstream Artifacts are linked to
        their real parent Artifact nodes with the proper Type information.
        '''
        all_metadata = dict()
        artifacts_as_metadata = []
        if (all_params := action_details.get('parameters')) is not None:
            for param in all_params:
                param_val, = param.values()
                if isinstance(param_val, MetadataInfo):
                    param_name, = param.keys()
                    md_fp = param_val.relative_fp
                    all_metadata.update({param_name: md_fp})

                    artifacts_as_metadata += [
                        {'artifact_passed_as_metadata': uuid} for uuid in
                        param_val.input_artifact_uuids
                    ]

        return all_metadata, artifacts_as_metadata

    def _parse_metadata(
        self, zf: ZipFile, metadata_fps: Dict[str, str]
    ) -> Dict[str, pd.DataFrame]:
        '''
        Parses all metadata files captured from Metadata and MetadataColumns
        (identifiable by !metadata tags) into pd.DataFrames.

        Parameters
        ----------
        zf : ZipFile
            The zipfile object of the archive.
        metadata_fps : dict
            A dict of parameter names to metadata filenames for metadata
            paramters.

        Returns
        -------
        dict
            A dict of parameter names to dataframe objects that is loaded from
            the corresponding metadata file.

            An empty dict if there is no metadata.
        '''
        if metadata_fps == {}:
            return {}

        root_uuid = get_root_uuid(zf)
        pfx = pathlib.Path(root_uuid) / 'provenance'
        if root_uuid == self._uuid:
            pfx = pfx / 'action'
        else:
            pfx = pfx / 'artifacts' / self._uuid / 'action'

        all_md = dict()
        for param_name in metadata_fps:
            filepath = str(pfx / metadata_fps[param_name])
            with zf.open(filepath) as fh:
                df = pd.read_csv(BytesIO(fh.read()), sep='\t')
                all_md[param_name] = df

        return all_md

    def __repr__(self) -> str:
        return repr(self._result_md)

    __str__ = __repr__

    def __hash__(self) -> int:
        return hash(self._uuid)

    def __eq__(self, other) -> bool:
        return (
            self.__class__ == other.__class__ and self._uuid == other._uuid
        )


class _Action:
    '''Provenance data from action.yaml for a single QIIME2 Result.'''

    @property
    def action_id(self) -> str:
        '''The UUID of the Action itself.'''
        return self._execution_details['uuid']

    @property
    def action_type(self) -> str:
        '''
        The type of Action represented e.g. Method, Pipeline, et al.
        '''
        return self._action_details['type']

    @property
    def runtime(self) -> timedelta:
        '''The elapsed run time of the Action, as a datetime object.'''
        end = self._execution_details['runtime']['end']
        start = self._execution_details['runtime']['start']
        return end - start

    @property
    def runtime_str(self) -> str:
        ''' The elapsed run time of the Action in seconds and microseconds.'''
        return self._execution_details['runtime']['duration']

    @property
    def action_name(self) -> str:
        '''
        The name of the action itself. Imports return 'import'.
        '''
        if self.action_type == 'import':
            return 'import'
        return self._action_details.get('action')

    @property
    def plugin(self) -> str:
        '''
        The plugin which executed this Action. Returns 'framework' if this is
        an import.
        '''
        if self.action_type == 'import':
            return 'framework'

        plugin = self._action_details.get('plugin')
        return plugin.replace('-', '_')

    @property
    def inputs(self) -> dict:
        '''
        Creates a dict of artifact inputs to this action.

        Returns
        -------
        dict
            A mapping of input name to uuid for each input passed to the
            corresponding action.

        Raises
        ------
        ValueError
            If an input ResultCollection is detected. These are detectable here
            because typically the inputs section is structured as:

            inputs:
            - some_input_name: some_uuid
            - some_other_input_name: some_other_uuid

            but when a ResultCollection is an input it is structed as:

            inputs:
            - result_collection_name:
                - some_key: some_uuid
                - some_other_key: some_other_uuid

            and thus is a different structure entirely.

        Notes
        -----
            When support for ResultCollections exists, the erroring should
            obviously be removed, and the way that items are added to `results`
            might need to be rethought--should nested dicts be added, or should
            ResultCollections be unpacked and and added piece by piece?
        '''
        inputs = self._action_details.get('inputs')
        results = {}
        if inputs is not None:
            for item in inputs:
                first_value = next(iter(item.values()))
                if type(first_value) is list:
                    msg = (
                        'An action in provenance took a ResultCollection '
                        'as input. Replay of ResultCollections is not '
                        'currently supported.'
                    )
                    raise ValueError(msg)

                results.update(item.items())

        return results

    @property
    def parameters(self) -> dict:
        '''Returns a dict of parameters passed to this action.'''
        params = self._action_details.get('parameters')
        results = {}
        if params is not None:
            for item in params:
                results.update(item.items())
        return results

    @property
    def output_name(self) -> Optional[str]:
        '''
        Gets the output name of the node.

        Returns
        -------
        str or None
            The name of the output as parsed from action.yaml, or None if there
            is no output-name section.

        Raises
        ------
        ValueError
            If the artifact comes from a ResultCollection. This is detectable
            here because outputs from a ResultCollection look like:

            output-name:
            - output
            - key
            - position/total positions

            and are thus a different structure entirely.

        Notes
        -----
            The erroring should be removed when support for ResultCollections
            are added.
        '''
        output_name = self._action_details.get('output-name')

        if type(output_name) is list:
            msg = (
                'A ResultCollection was returned by an action in provenance. '
                'Replay of ResultCollections are not currently supported.'
            )
            raise ValueError(msg)

        return output_name

    @property
    def format(self) -> Optional[str]:
        '''Returns this action's format field if any.'''
        return self._action_details.get('format')

    @property
    def transformers(self) -> Optional[Dict]:
        '''Returns this action's transformers dictionary if any.'''
        return self._action_dict.get('transformers')

    def __init__(self, zf: ZipFile, fp: str):
        with tempfile.TemporaryDirectory() as tempdir:
            zf.extractall(tempdir)
            action_fp = os.path.join(tempdir, fp)
            with open(action_fp) as fh:
                self._action_dict = yaml.safe_load(fh)

        self._action_details = self._action_dict['action']
        self._execution_details = self._action_dict['execution']

    def __repr__(self):
        return (
            f'_Action(action_id={self.action_id}, type={self.action_type},'
            f' plugin={self.plugin}, action={self.action_name})'
        )


class _Citations:
    '''
    Citations for a single QIIME2 Result, as a dict of citation dicts keyed
    on the citation's bibtex ID.
    '''
    def __init__(self, zf: ZipFile, fp: str):
        bib_db = bp.loads(zf.read(fp))
        self.citations = bib_db.get_entry_dict()

    def __repr__(self):
        keys = list(self.citations.keys())
        return f'Citations({keys})'


class _ResultMetadata:
    '''Basic metadata about a single QIIME2 Result from metadata.yaml.'''
    def __init__(self, zf: ZipFile, md_fp: str):
        _md_dict = yaml.safe_load(zf.read(md_fp))
        self.uuid = _md_dict['uuid']
        self.type = _md_dict['type']
        self.format = _md_dict['format']

    def __repr__(self):
        return (
            f'UUID:\t\t{self.uuid}\n'
            f'Type:\t\t{self.type}\n'
            f'Data Format:\t{self.format}'
        )


class Parser(metaclass=abc.ABCMeta):
    accepted_data_types: str

    @classmethod
    @abc.abstractmethod
    def get_parser(cls, artifact_data: Any) -> 'Parser':
        '''
        Return the appropriate Parser if this Parser type can handle the data
        passed in.

        Should raise an appropriate exception if this Parser cannot handle the
        data.
        '''

    @abc.abstractmethod
    def parse_prov(self, cfg: Config, data: Any) -> ParserResults:
        '''
        Parse provenance to return a ParserResults.
        '''


class ArchiveParser(Parser):
    accepted_data_types = 'a path to a file (a string) or a file-like object'

    @classmethod
    def get_parser(cls, artifact: Union[str, pathlib.PosixPath]) -> Parser:
        '''
        Returns the correct archive format parser for a zip archive.

        Parameters
        ----------
        artifact_data : str or pathlib.PosixPath
            A path to a zipped archive.

        Returns
        -------
        Parser
            An ArchiveParser object for the version of the artifact. One of
            ParserV[0-6].
        '''
        if type(artifact) is pathlib.PosixPath:
            artifact = str(artifact)
        if type(artifact) is not str:
            raise TypeError(
                'ArchiveParser expects a string or pathlib.PosixPath path to '
                f'an archive, not an object of type {str(type(artifact))}.'
            )
        if os.path.isdir(artifact):
            raise ValueError('ArchiveParser expects a file, not a directory.')

        try:
            with ZipFile(artifact, 'r') as zf:
                archive_version, _ = parse_version(zf)
            return FORMAT_REGISTRY[archive_version]()
        except KeyError as e:
            raise KeyError(
                f'While trying to parse artifact {artifact}, '
                'a corresponding parser was not found for archive version '
                f'{archive_version}: {str(e)}.'
            )

    def parse_prov(cls, cfg: Config, data: Any) -> ParserResults:
        raise NotImplementedError(
            'Use a subclass that usefully defines parse_prov for some format.'
        )


class ParserV0(ArchiveParser):
    '''
    Parser for V0 archives. V0 archives have no ancestral provenance.
    '''
    # These are files we expect will be present in every QIIME2 archive with
    # this format. "Optional" filenames (like Metadata, which may or may
    # not be present in an archive) should not be included here.
    expected_files_root_only = tuple()
    expected_files_all_nodes = ('metadata.yaml', 'VERSION')

    def parse_prov(self, cfg: Config, archive: str) -> ParserResults:
        '''
        Parses an artifact's provenance into a directed acyclic graph.

        In the case of v0 archives, the only provenance information is that
        which is attached to the artifact itself; information about ancestor
        nodes does not exist. The parsed dag contains only a single node.

        In the case of v1 archives, ancestor nodes do exist in the
        archive. However, because the corresponding action.yaml does not track
        output names, when two outputs share the same semantic type, it is not
        possible to untangle provenance. Instead of wrangling with this and
        in consideration of the expected rarity of v1 archives, it was decided
        to treat v1 archives as v0 archives.

        Parameters
        ----------
        cfg : Config
            A dataclass that stores four boolean flags: whether to perform
            checksum validation, whether to parse study metadata, whether to
            recursively parse nested directories, and whether to enable verbose
            mode.
        archive : str
            A path to the artifact to be parsed.

        Returns
        -------
        ParserResults
            A dataclass that stores the parsed artifact uuids, the parsed
            networkx graph, the provenance-is-valid flag, and the
            checksum diff.
        '''
        with ZipFile(archive) as zf:
            if cfg.perform_checksum_validation:
                provenance_is_valid, checksum_diff = \
                    self._validate_checksums(zf)
            else:
                provenance_is_valid = ValidationCode.VALIDATION_OPTOUT
                checksum_diff = None

            root_uuid = get_root_uuid(zf)
            warnings.warn(
                f'Artifact {root_uuid} was created prior to provenance '
                'tracking. Provenance data will be incomplete.',
                UserWarning
            )

            exp_node_fps = []
            for fp in self.expected_files_all_nodes:
                exp_node_fps.append(pathlib.Path(root_uuid) / fp)

            prov_fps = self._get_provenance_fps(zf)
            self._assert_expected_files_present(
                zf, exp_node_fps, prov_fps
            )
            # we have confirmed that all expected fps for this node exist
            node_fps = exp_node_fps

            nodes = {}
            nodes[root_uuid] = ProvNode(cfg, zf, node_fps)
            graph = self._digraph_from_archive_contents(nodes)

        return ParserResults(
            {root_uuid},
            graph,
            provenance_is_valid,
            checksum_diff
        )

    def _parse_root_md(self, zf: ZipFile, root_uuid: str) -> _ResultMetadata:
        '''
        Parses the root metadata file of an archive for its uuid, semantic
        type, and format.

        Parameters
        ----------
        zf : ZipFile
            A zipfile object of a v0 artifact.
        root_uuid : str
            The uuid of the root node. Because this operates on a v0 archive,
            the root node is the only node.

        Returns
        -------
        _ResultMetadata
            An object representing the information stored in a metadata.yaml
            file, namely the uuid, type, and format fields.
        '''
        root_md_fp = os.path.join(root_uuid, 'metadata.yaml')
        if root_md_fp not in zf.namelist():
            raise ValueError(
                'Malformed Archive: root metadata.yaml file '
                f'misplaced or nonexistent in {zf.filename}'
            )
        return _ResultMetadata(zf, root_md_fp)

    def _validate_checksums(
            self, zf: ZipFile
    ) -> Tuple[ValidationCode, Optional[ChecksumDiff]]:
        '''
        Return the ValidationCode and ChecksumDiff for an archive. Because
        checksums were not introduced, until ArchiveFormat version 5,
        uses the PREDATES_CHECKSUMS flag and returns None to indicate that
        checksum diffing was not performed.

        Parameters
        ----------
        zf : ZipFile
            The zipfile object representing the archive. Ignored here but
            needed in signature for inheritance.

        Returns
        -------
        tuple of (ValidationCode, None)
            The validation code and None to indicate missing ChecksumDiff.
        '''
        return (ValidationCode.PREDATES_CHECKSUMS, None)

    def _digraph_from_archive_contents(
        self, archive_contents: Dict[str, 'ProvNode']
    ) -> nx.DiGraph:
        '''
        Builds a networkx.DiGraph from a {UUID: ProvNode} dictionary.

        1. Create an empty nx.digraph.
        2. Gather nodes and their required attributes and add them to the
           DiGraph.
        3. Add edges to graph (including all !no-provenance nodes)
        4. Create guaranteed node attributes for these no-provenance nodes,
           which wouldn't otherwise have them.

        Parameters
        ----------
        archive_contents : dict of {str to ProvNode}
            A dictionary of node uuids to their representative ProvNode
            objects.

        Returns
        -------
        nx.DiGraph
            The directed, acyclic graph representation of the provenance of
            the archive. Edge directionality is from parent to child. Parents
            may have multiple children and children may have multiple parents.
        '''
        dag = nx.DiGraph()
        nodes = []
        for node_uuid, node in archive_contents.items():
            node_info = {
                'node_data': node,
                'has_provenance': node.has_provenance
            }
            nodes.append((node_uuid,  node_info))
        dag.add_nodes_from(nodes)

        edges = []
        for node_uuid, attrs in dag.nodes(data=True):
            if parents := attrs['node_data']._parents:
                for parent in parents:
                    parent_uuid, = parent.values()
                    edges.append((parent_uuid, node_uuid))
        dag.add_edges_from(edges)

        return dag

    def _get_provenance_fps(self, zf: ZipFile) -> List[pathlib.Path]:
        '''
        Collect filepaths of all provenance-relevant files in an archive.
        Relevant is defined by `self.expected_files_all_nodes` and
        `self.expected_files_root_only` (which is empty).

        Parameters
        ----------
        zf : ZipFile
            The zipfile object of the archive.

        Returns
        -------
        list of pathlib.Path
            Filepaths relative to root of zipfile for each file of interest.
        '''
        fps = []
        for fp in zf.namelist():
            for expected_filename in self.expected_files_all_nodes:
                if expected_filename in fp:
                    fps.append(pathlib.Path(fp))

        return fps

    def _assert_expected_files_present(
            self,
            zf: ZipFile,
            expected_node_fps: List[pathlib.Path],
            prov_fps: List[pathlib.Path],
    ):
        '''
        Makes sure that all expected files for a given node are present in an
        archive. Raises a ValueError if not.

        Parameters
        ----------
        zf : ZipFile
            The zipfile object representing an archive.
        expected_node_fps : list of pathlib.Path
            The filepaths that are expected to be present in the zipfile for
            some node.
        prov_fps : list of pathlib.Path
            All provenance-relevant filepaths in the archive.

        Raises
        ------
        ValueError
            If there are expected provenance-relevant files missing from
            a node.
        '''
        error_contents = 'Malformed Archive: '
        root_uuid = get_root_uuid(zf)
        for fp in expected_node_fps:
            if fp not in prov_fps:
                node_uuid = get_nonroot_uuid(fp)
                error_contents += (
                    f'{fp.name} file for node {node_uuid} '
                    f'misplaced or nonexistent in {zf.filename}.\n'
                )
                error_contents += (
                    f'Archive {root_uuid} may be corrupt '
                    'or provenance may be false.'
                )
                raise ValueError(error_contents)


class ParserV1(ParserV0):
    '''
    Parser for V1 archives. Although action.yaml was introduced for this
    archive version, we are pretending that it was introduced in V2 because of
    difficulties untangling provenance without output names. V1 archives are
    treated as having no provenance, like V0 archives.
    '''
    expected_files_root_only = ParserV0.expected_files_root_only
    expected_files_all_nodes = ParserV0.expected_files_all_nodes


class ParserV2(ParserV1):
    '''
    Parser for V2 archives. Introduces action/action.yaml to provenance.
    Directory structure identical to V1, action.yaml changes to support
    Pipelines.
    '''
    expected_files_root_only = ParserV1.expected_files_root_only
    expected_files_all_nodes = (
        *ParserV1.expected_files_all_nodes, 'action/action.yaml'
    )

    def parse_prov(self, cfg: Config, archive: str) -> ParserResults:
        '''
        Parses an artifact's provenance into a directed acyclic graph.

        For each artifact in provenance, gathers all corresponding
        provenance-relevant files and constructs a ProvNode. Once all
        ProvNodes are constructed, creates the provenance graph.

        Parameters
        ----------
        cfg : Config
            A dataclass that stores four boolean flags: whether to perform
            checksum validation, whether to parse study metadata, whether to
            recursively parse nested directories, and whether to enable verbose
            mode.
        archive_data : str
            A path to the artifact to be parsed.

        Returns
        -------
        ParserResults
            A dataclass that stores the parsed artifact uuids, the parsed
            networkx graph, the provenance-is-valid flag, and the
            checksum diff.
        '''
        with ZipFile(archive) as zf:
            if cfg.perform_checksum_validation:
                provenance_is_valid, checksum_diff = \
                    self._validate_checksums(zf)
            else:
                provenance_is_valid = ValidationCode.VALIDATION_OPTOUT
                checksum_diff = None

            prov_fps = self._get_provenance_fps(zf)
            root_uuid = get_root_uuid(zf)

            # make a provnode for each UUID
            archive_contents = {}
            for fp in prov_fps:
                exp_node_fps = []
                if 'artifacts' not in fp.parts:
                    node_uuid = root_uuid
                    prefix = pathlib.Path(node_uuid) / 'provenance'
                    root_only_expected_fps = []
                    for exp_filename in self.expected_files_root_only:
                        root_only_expected_fps.append(
                            pathlib.Path(node_uuid) / exp_filename
                        )
                    exp_node_fps += root_only_expected_fps
                else:
                    node_uuid = get_nonroot_uuid(fp)
                    # /root-uuid/provenance/artifacts/node-uuid
                    prefix = pathlib.Path(*fp.parts[0:4])

                if node_uuid in archive_contents:
                    continue

                for expected_file in self.expected_files_all_nodes:
                    exp_node_fps.append(prefix / expected_file)

                self._assert_expected_files_present(
                    zf, exp_node_fps, prov_fps
                )
                # we have confirmed that all expected fps for this node exist
                node_fps = exp_node_fps

                archive_contents[node_uuid] = ProvNode(cfg, zf, node_fps)

        graph = self._digraph_from_archive_contents(archive_contents)

        return ParserResults(
            {root_uuid},
            graph,
            provenance_is_valid,
            checksum_diff
        )

    def _get_provenance_fps(self, zf: ZipFile) -> List[pathlib.Path]:
        '''
        Collect filepaths of all provenance-relevant files in an archive.
        Relevant is defined by `self.expected_files_all_nodes` and
        `self.expected_files_root_only`.

        Parameters
        ----------
        zf : ZipFile
            The zipfile object of the archive.

        Returns
        -------
        list of pathlib.Path
            Filepaths relative to root of zipfile for each file of interest.
        '''
        fps = []
        for fp in zf.namelist():
            for expected_filename in self.expected_files_all_nodes:
                if 'provenance' in fp and expected_filename in fp:
                    fps.append(pathlib.Path(fp))

        root_uuid = get_root_uuid(zf)
        for expected_filename in self.expected_files_root_only:
            fps.append(pathlib.Path(root_uuid) / expected_filename)

        return fps


class ParserV3(ParserV2):
    '''
    Parser for V3 archives. Directory structure identical to V1 & V2,
    action.yaml now supports variadic inputs, so !set tags in action.yaml.
    '''
    expected_files_root_only = ParserV2.expected_files_root_only
    expected_files_all_nodes = ParserV2.expected_files_all_nodes


class ParserV4(ParserV3):
    '''
    Parser for V4 archives. Adds citations to directory structure, changes to
    action.yaml including transformers.
    '''
    expected_files_root_only = ParserV3.expected_files_root_only
    expected_files_all_nodes = (
        *ParserV3.expected_files_all_nodes, 'citations.bib'
    )


class ParserV5(ParserV4):
    '''
    Parser for V5 archives. Adds checksum validation with checksums.md5.
    '''
    expected_files_root_only = ('checksums.md5', )
    expected_files_all_nodes = ParserV4.expected_files_all_nodes

    def _validate_checksums(
            self, zf: ZipFile
    ) -> Tuple[ValidationCode, Optional[ChecksumDiff]]:
        '''
        Checksum support added for v5, so perform checksum validation.

        Parameters
        ----------
        zf : ZipFile
            The zipfile object representation of the parsed archive.

        Returns
        -------
        tuple of (ValidationCode, ChecksumDiff or None)
            Where ValidationCode is one of valid, invalid, predates checksums,
            optout.
            Where ChecksumDiff contains filepaths of all changed, added, and
            removed files since last checksumming.
            If checksums.md5 is missing from archive the archive, an invalid
            code is returned and a ChecksumDiff of None is returned.

        Notes
        -----
        Because a ChecksumDiff of None here has a different interpetation
        than in pre-V5 archive parsers, the ChecksumDiff should only be
        intepreted in conjuction with the ValidationCode.
        '''
        return validate_checksums(zf)


class ParserV6(ParserV5):
    '''
    Parser for V6 archives. Adds support for output collections, adds
    execution_context field to action.yaml.
    '''
    expected_files_root_only = ParserV5.expected_files_root_only
    expected_files_all_nodes = ParserV5.expected_files_all_nodes


FORMAT_REGISTRY = {
    # NOTE: update for new format versions in qiime2.core.archive.Archiver
    '0': ParserV0,
    '1': ParserV1,
    '2': ParserV2,
    '3': ParserV3,
    '4': ParserV4,
    '5': ParserV5,
    '6': ParserV6
}
