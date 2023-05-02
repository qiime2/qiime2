from typing import Any, List, NamedTuple, Set, Union
import warnings

from .util import UUID


def citation_key_constructor(loader, node) -> str:
    """
    A constructor for !cite yaml tags, returning a bibtex key as a str.
    All we need for now is a key string we can match in citations.bib,
    so _we're not parsing these into component substrings_.

    If that need arises in future, these are spec'ed in provenance.py as:
    <domain>|<package>:<version>|[<identifier>|]<index>

    and frequently look like this (note no identifier):
    framework|qiime2:2020.6.0.dev0|0
    """
    value = loader.construct_scalar(node)
    return value


def color_constructor(loader, node) -> str:
    """
    Constructor for !color tags, returning an str.
    Color was a primitive type representing a 3 or 6 digit color hex code,
    matching ^#(?:[0-9a-fA-F]{3}){1,2}$

    Per E. Bolyen,these were unused by any plugins. They were removed in
    e58ed5f8ba453035169d560e0223e6a37774ae08, which was released in 2019.4
    """
    return loader.construct_scalar(node)


class MetadataInfo(NamedTuple):
    """
    A namedtuple representation of the data in one !metadata yaml tag.
    """
    input_artifact_uuids: List[UUID]
    relative_fp: str


def metadata_path_constructor(loader, node) -> MetadataInfo:
    """
    A constructor for !metadata yaml tags, which come in the form
    [<uuid_ref>[,<uuid_ref>][...]:]<relative_filepath>

    Returns a MetadataInfo object containing a list of UUIDs, and the relative
    filepath where the metadata was written into the zip archive

    Most commonly, we see:
    !metadata 'sample_metadata.tsv'

    In cases where Artifacts are used as metadata, we see:
    !metadata '415409a4-371d-4c69-9433-e3eaba5301b4:feature_metadata.tsv'

    In cases where multiple Artifacts as metadata were merged,
    it is possible for multiple comma-separated uuids to precede the ':'
    !metadata '<uuid1>,<uuid2>,...,<uuidn>:feature_metadata.tsv'

    The metadata files (including "Artifact metadata") are saved in the same
    dir as `action.yaml`. The UUIDs listed must be incorporated into our
    provenance graph as parents, so are returned in list form.
    """
    raw = loader.construct_scalar(node)
    if ':' in raw:
        artifact_uuids, rel_fp = raw.split(':')
        artifact_uuids = artifact_uuids.split(',')
    else:
        artifact_uuids = []
        rel_fp = raw
    return MetadataInfo(artifact_uuids, rel_fp)


def no_provenance_constructor(loader, node) -> UUID:
    """
    Constructor for !no-provenance tags. These tags are written by QIIME 2 when
    an input has no /provenance dir, as in the case of v0 archives that have
    been used in analyses in QIIME2 V1+. They look like this:

    action:
       inputs:
       -   table: !no-provenance '34b07e56-27a5-4f03-ae57-ff427b50aaa1'

    For now at least, this constructor warns but otherwise disregards the
    no-provenance-ness of these. The v0 parser deals with them directly anyway.
    """
    uuid = loader.construct_scalar(node)
    warnings.warn(f"Artifact {uuid} was created prior to provenance tracking. "
                  + "Provenance data will be incomplete.", UserWarning)
    return uuid


def ref_constructor(loader, node) -> Union[str, List[str]]:
    """
    A constructor for !ref yaml tags. These tags describe yaml values that
    reference other namespaces within the document, using colons to separate
    namespaces. For example:
    !ref 'environment:plugins:sample-classifier'

    At present, ForwardRef tags are only used in the framework to 'link' the
    plugin name to the plugin version and other details in the 'execution'
    namespace of action.yaml

    This constructor explicitly handles this type of !ref by extracting and
    returning the plugin name to simplify parsing, while supporting the return
    of a generic list of 'keys' (e.g. ['environment', 'framework', 'version'])
    in the event ForwardRef is used more broadly in future.
    """
    value = loader.construct_scalar(node)
    keys = value.split(':')
    if keys[0:2] == ['environment', 'plugins']:
        plugin_name = keys[2]
        return plugin_name
    else:
        return keys


def set_constructor(loader, node) -> Set[Any]:
    """
    A constructor for !set yaml tags, returning a python set object
    """
    value = loader.construct_sequence(node)
    return set(value)


# NOTE: New yaml tag constructors must be added to this registry, or tags will
# raise ConstructorErrors
CONSTRUCTOR_REGISTRY = {
    '!cite': citation_key_constructor,
    '!color': color_constructor,
    '!metadata': metadata_path_constructor,
    '!no-provenance': no_provenance_constructor,
    '!ref': ref_constructor,
    '!set': set_constructor,
    }
