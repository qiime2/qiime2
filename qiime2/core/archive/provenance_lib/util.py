import pathlib
import re
import zipfile

# Alias string as UUID so we can specify types more clearly
UUID = str

# Alias string as FileName because that's what some strings mean
# FileNames are not path objects - just strings that describe paths
FileName = str

### NOTE: implemented as method on _Archive class
def get_root_uuid(zf: zipfile.ZipFile) -> UUID:
    """
    Returns the root UUID for a QIIME 2 Archive.

    There's no particular reason we use the first filename here.  All QIIME 2
    Artifacts store their contents in a directory named with the Artifact's
    UUID, so we can get the UUID of an artifact by taking the first part of the
    filepath of any file in the zip archive.
    """
    return pathlib.Path(zf.namelist()[0]).parts[0]


def get_nonroot_uuid(fp: pathlib.Path) -> UUID:
    """
    For non-root provenance files, get the Result's uuid from the path
    (avoiding the root Result's UUID which is in all paths)
    """
    if fp.name == 'action.yaml':
        uuid = fp.parts[-3]
    else:
        uuid = fp.parts[-2]
    return uuid


def camel_to_snake(name: str) -> str:
    """
    There are more comprehensive and faster ways of doing this (incl compiling)
    but it handles acronyms in semantic types nicely
    e.g. EMPSingleEndSequences -> emp_single_end_sequences
    c/o https://stackoverflow.com/a/1176023/9872253
    """
    # this will frequently be called on QIIME type expressions, so drop [ and ]
    name = re.sub(r'[\[\]]', '', name)
    # camel to snake
    name = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', name).lower()
