# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import contextlib
import warnings
import hashlib
import stat
import os
import io
import re
import collections
import uuid as _uuid
import yaml
import zipfile
import pathlib
import shutil
import subprocess

import decorator

from qiime2.core.annotate import ANNOTATION_TYPE_LIST, UnknownAnnotation

READ_ONLY_FILE = stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH
READ_ONLY_DIR = stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | stat.S_IRUSR \
    | stat.S_IRGRP | stat.S_IROTH
USER_GROUP_RWX = stat.S_IRWXU | stat.S_IRWXG
OTHER_NO_WRITE = stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH
# public key algorithm identifiers
# https://datatracker.ietf.org/doc/html/rfc4880#section-9.1
_PUBKEY_ALG = {
    '1': 'RSA',
    '2': 'RSA',
    '3': 'RSA',
    '16': 'ElGamal',
    '17': 'DSA',
    '18': 'ECDH',
    '19': 'ECDSA',
    '22': 'EdDSA'
}


def get_view_name(view):
    from .format import FormatBase
    if not isinstance(view, type):
        view = view.__class__

    if issubclass(view, FormatBase):
        # Not qualname because we don't have a notion of "nested" formats
        return view.__name__

    return ':'.join([view.__module__, view.__qualname__])


def tuplize(x):
    if type(x) is not tuple:
        return (x,)
    return x


def overrides(cls):
    def decorator(func):
        if not hasattr(cls, func.__name__):
            raise AssertionError("%r does not override %r"
                                 % (func, cls.__name__))
        return func
    return decorator


def superscript(number):
    table = {
        '0': chr(8304), '1': chr(185), '2': chr(178), '3': chr(179),
        **{str(i): chr(x) for i, x in enumerate(range(8308, 8314), 4)},
        'a': chr(7491), 'e': chr(7497), 'f': chr(7584), 'i': chr(8305),
        'n': chr(8319), '-': chr(8315), '.': chr(39), ',': chr(39)
    }
    return ''.join([table[d] for d in str(number)])


def find_duplicates(iterable):
    """Find duplicate values in an iterable.

    Parameters
    ----------
    iterable : iterable
        Iterable to search for duplicates.

    Returns
    -------
    set
        Values that are duplicated in `iterable`.

    Notes
    -----
    Values in `iterable` must be hashable.

    """
    # Modified from https://stackoverflow.com/a/9835819/3776794 to return
    # duplicates instead of remove duplicates from an iterable.
    seen = set()
    duplicates = set()
    for value in iterable:
        if value in seen:
            duplicates.add(value)
        else:
            seen.add(value)
    return duplicates


# Concept from: http://stackoverflow.com/a/11157649/579416
def duration_time(relative_delta):
    attrs = ['years', 'months', 'days', 'hours', 'minutes', 'seconds',
             'microseconds']
    results = []
    for attr in attrs:
        value = getattr(relative_delta, attr)
        if value != 0:
            if value == 1:
                # Remove plural 's'
                attr = attr[:-1]
            results.append("%d %s" % (value, attr))
    if results:
        text = results[-1]
        if results[:-1]:
            text = ', and '.join([', '.join(results[:-1]), text])
        return text
    else:
        # Great Scott! No time has passed!
        return '0 %s' % attrs[-1]


def has_checksum_native(checksum_type):
    if checksum_type == 'md5':
        checksum_util = 'md5sum'
    elif checksum_type == 'sha512':
        checksum_util = 'sha512'

    return shutil.which(f'{checksum_util}') is not None


def checksum(filepath, checksum_type):
    if has_checksum_native(checksum_type):
        return checksum_native(filepath, checksum_type)
    else:
        return checksum_python(filepath, checksum_type)


def checksum_python(filepath, checksum_type):
    if checksum_type == 'md5':
        hash_obj = hashlib.md5()
    elif checksum_type == 'sha512':
        hash_obj = hashlib.sha512()
    # we shouldn't ever hit this branch, but just in case
    else:
        raise TypeError(f'Unsupported checksum type: {checksum_type!r}')

    with open(str(filepath), mode='rb') as fh:
        for chunk in iter(lambda: fh.read(io.DEFAULT_BUFFER_SIZE), b""):
            hash_obj.update(chunk)
    return hash_obj.hexdigest()


def checksum_native(filepath, checksum_type):
    if checksum_type == 'md5':
        cmd = ['md5sum', str(filepath)]
    elif checksum_type == 'sha512':
        cmd = ['sha512sum', str(filepath)]
    # we shouldn't ever hit this branch, but just in case
    else:
        raise TypeError(f'Unsupported checksum type: {checksum_type!r}')

    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    _, digest = from_checksum_format(result.stdout)
    return digest


def checksum_zip(zf: zipfile.ZipFile, filepath: str,
                 checksum_type: str) -> str:
    """
    Given a ZipFile object and relative filepath within the zip archive,
    returns the checksum of the file
    """
    if checksum_type == 'md5':
        hash_obj = hashlib.md5()
    elif checksum_type == 'sha512':
        hash_obj = hashlib.sha512()
    # we shouldn't ever hit this branch, but just in case
    else:
        raise TypeError(f'Unsupported checksum type: {checksum_type!r}')

    with zf.open(filepath) as fh:
        for chunk in iter(lambda: fh.read(io.DEFAULT_BUFFER_SIZE), b""):
            hash_obj.update(chunk)
    return hash_obj.hexdigest()


def checksum_directory(directory, checksum_type):
    if has_checksum_native(checksum_type):
        checksum = checksum_native
    else:
        checksum = checksum_python

    directory = str(directory)
    sums = collections.OrderedDict()
    for root, dirs, files in os.walk(directory, topdown=True):
        dirs[:] = sorted([d for d in dirs if not d[0] == '.'])
        for file in sorted(files):
            if file[0] == '.':
                continue

            path = os.path.join(root, file)
            sums[os.path.relpath(path, start=directory)] = \
                checksum(path, checksum_type)
    return sums


def checksum_directory_zip(zf: zipfile.ZipFile, checksum_type: str) -> dict:
    """
    Returns a mapping of fp/checksum pairs for all files in zf.

    The root dir has been removed from these filepaths. This mimics the output
    in checksums.md5 (without sorted descent), but is not generalizable beyond
    QIIME 2 archives.
    """
    sums = dict()
    for file in zf.namelist():
        fp = pathlib.Path(file)
        if fp.name == f'checksums.{checksum_type}':
            continue

        file_parts = list(fp.parts)
        internal_path_parts = file_parts[1:]

        if internal_path_parts and internal_path_parts[0] == 'annotations':
            continue

        fp_w_o_root_uuid = pathlib.Path(*(file_parts[1:]))
        sums[str(fp_w_o_root_uuid)] = checksum_zip(zf, file, checksum_type)

    return sums


def to_checksum_format(filepath, checksum):
    # see https://www.gnu.org
    # /software/coreutils/manual/html_node/md5sum-invocation.html
    if '\\' in filepath or '\n' in filepath:
        filepath = filepath.replace('\\', '\\\\').replace('\n', '\\n')
        checksum = '\\' + checksum

    return '%s  %s' % (checksum, filepath)


def from_checksum_format(line):
    line = line.rstrip('\n')
    parts = line.split('  ', 1)
    if len(parts) < 2:
        parts = line.split(' *', 1)

    checksum, filepath = parts

    if checksum[0] == '\\':
        chars = ''
        escape = False
        # Gross, but regular `.replace` will overlap with itself and
        # negative lookbehind in regex is *probably* harder than scanning
        for char in filepath:
            # 1) Escape next character
            if not escape and char == '\\':
                escape = True
                continue

            # 2) Handle escape sequence
            if escape:
                try:
                    chars += {'\\': '\\', 'n': '\n'}[char]
                except KeyError:
                    chars += '\\' + char  # Wasn't an escape after all
                escape = False
                continue

            # 3) Nothing interesting
            chars += char

        checksum = checksum[1:]
        filepath = chars

    return filepath, checksum


@contextlib.contextmanager
def warning():
    def _warnformat(msg, category, filename, lineno, file=None, line=None):
        return '%s:%s: %s: %s\n' % (filename, lineno, category.__name__, msg)

    default_warn_format = warnings.formatwarning
    try:
        warnings.formatwarning = _warnformat
        warnings.filterwarnings('always')
        yield warnings.warn
    finally:
        warnings.formatwarning = default_warn_format


# Descriptor protocol for creating an attribute that is bound to an
# (arbitrarily nested) attribute accessible to the instance at runtime.
class LateBindingAttribute:
    def __init__(self, attribute):
        self._attribute = attribute

    def __get__(self, obj, cls=None):
        attrs = self._attribute.split('.')
        curr_attr = obj
        for attr in attrs:
            curr_attr = getattr(curr_attr, attr)
        return staticmethod(curr_attr).__get__(obj, cls)


# Removes the first parameter from a callable's signature.
class DropFirstParameter(decorator.FunctionMaker):
    @classmethod
    def from_function(cls, function):
        return cls.create(function, "return None", {})

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.signature = self._remove_first_arg(self.signature)
        self.shortsignature = self._remove_first_arg(self.shortsignature)

    def _remove_first_arg(self, string):
        return ",".join(string.split(',')[1:])[1:]


def _immutable_error(obj, *args):
    raise TypeError('%s is immutable.' % obj.__class__.__name__)


class ImmutableBase:
    def _freeze_(self):
        """Disables __setattr__ when called. It is idempotent."""
        self._frozen = True  # The particular value doesn't matter

    __delattr__ = __setitem__ = __delitem__ = _immutable_error

    def __setattr__(self, *args):
        # This doesn't stop silly things like
        # object.__setattr__(obj, ...), but that's a pretty rude thing
        # to do anyways. We are just trying to avoid accidental mutation.
        if hasattr(self, '_frozen'):
            _immutable_error(self)
        super().__setattr__(*args)


def sorted_poset(iterable, *, key=None, reverse=False):
    values = list(iterable)
    elements = values
    if key is not None:
        elements = [key(x) for x in values]

    result = []
    sorted_elements = []
    for value, element in zip(values, elements):
        idx = 0
        for idx, placed in enumerate(sorted_elements, 1):
            if element <= placed:
                idx -= 1
                break

        result.insert(idx, value)
        sorted_elements.insert(idx, element)
    if reverse:
        result = list(reversed(result))
    return result


def is_uuid4(uuid_str):
    # Adapted from https://gist.github.com/ShawnMilo/7777304
    try:
        uuid = _uuid.UUID(hex=uuid_str, version=4)
    except ValueError:
        # The string is not a valid hex code for a UUID.
        return False

    # If uuid_str is a valid hex code, but an invalid uuid4, UUID.__init__
    # will convert it to a valid uuid4.
    return str(uuid) == uuid_str


def set_permissions(path, file_permissions=None, dir_permissions=None,
                    skip_root=False):
    """Set permissions on all directories and files under and including path
    """
    # Panfs is currently causing issues for us setting permissions. We still
    # want to set rwx for user and group before we remove things to ensure we
    # can remove them, but we want to temporarily no-op other permission
    # changes
    if file_permissions != USER_GROUP_RWX:
        file_permissions = None

    if dir_permissions != USER_GROUP_RWX:
        dir_permissions = None

    # Just get out if we aren't doing anything
    if file_permissions is None and dir_permissions is None:
        return

    for directory, _, files in os.walk(path):
        # We may want to set permissions under a directory but not on the
        # directory itself.
        if dir_permissions and not (skip_root and directory == str(path)):
            try:
                os.chmod(directory, dir_permissions)
            except FileNotFoundError:
                pass

        for file in files:
            if file_permissions:
                try:
                    os.chmod(os.path.join(directory, file), file_permissions)
                except FileNotFoundError:
                    pass


def touch_under_path(path):
    """Touches everything under a given path to ensure they don't get culled by
    Mac
    """
    for directory, _, files in os.walk(path):
        try:
            os.utime(directory, None, follow_symlinks=False)
        except FileNotFoundError:
            pass

        for file in files:
            try:
                os.utime(
                    os.path.join(directory, file), None, follow_symlinks=False)
            except FileNotFoundError:
                pass


def load_action_yaml(path):
    """Takes a path to an unzipped Artifact and loads its action.yaml with
    yaml.safe_load
    """
    # TODO: Make these actually do something useful at least for the tags
    # that are relevant to what we need out of provenance (this is partially
    # done)
    def ref_constructor(loader, node):
        # We only care about the name of the thing we are referencing which
        # is at the end of this list
        return node.value.split(':')[-1]

    def cite_constructor(loader, node):
        return node.value

    def metadata_constructor(loader, node):
        # Use the checksum of the metadata as its identifier, so we can tell
        # if two artifacts used the same metadata input
        metadata_path = prov_path / node.value
        return checksum(filepath=metadata_path, checksum_type='md5')

    # these are backstops and are generally superceded by yaml.SafeLoader
    # which has the preferred constructors from provenance
    # found under CONSTRUCTOR_REGISTRY within provenance.py
    yaml.constructor.SafeConstructor.add_constructor('!ref', ref_constructor)
    yaml.constructor.SafeConstructor.add_constructor('!cite', cite_constructor)
    yaml.constructor.SafeConstructor.add_constructor('!metadata',
                                                     metadata_constructor)

    prov_path = path / 'provenance' / 'action'
    action_path = prov_path / 'action.yaml'

    with open(action_path) as fh:
        prov = yaml.safe_load(fh)

    return prov


def create_collection_name(*, name, key, idx, size):
    """ Only accepts kwargs. Creates a name for a collection item in a
        standardized way. Assumes 0 based indexing.
    """
    return [name, key, f'{idx + 1}/{size}']


# annotation helpers
# helper for parsing user name/email for keypair identification
def _parse_uid(uid_str):
    email_regex = re.compile(r'.*<([^>]+)>')
    uid_match = email_regex.match(uid_str or "")
    if uid_match:
        email = uid_match.group(1)
        name = uid_str[: uid_str.index('<')].strip()
        return name or None, email or None
    return (uid_str.strip() or None, None)


# helper for normalizing fingerprint formatting
def _normalize_fingerprint(s):
    return re.sub(r'\s+', '', (s or '')).upper()


def add_shared_attrs(annotation_type, meta_yaml):
    if annotation_type in ANNOTATION_TYPE_LIST:
        annotation = annotation_type.__new__(annotation_type)
    else:
        annotation = UnknownAnnotation.__new__(UnknownAnnotation)

    annotation.id = meta_yaml['id']
    annotation.name = meta_yaml['name']
    annotation.annotation_type = meta_yaml['type']
    annotation.created_at = meta_yaml['created_at']

    return annotation


# helper for locating root_fp for a given Result
def find_root_fp(annotations_dir, root_result_uuid):
    split_fp = annotations_dir.split(os.sep)
    root_result_uuid_index = split_fp.index(root_result_uuid)
    root_fp = os.sep.join(split_fp[0: root_result_uuid_index + 1])
    return root_fp


# helper for calculating the root level checksum digest
def sha512_file_hex(path):
    hex = hashlib.sha512()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1024*1024), b""):
            hex.update(chunk)
    return hex.hexdigest()


# helper for pulling keypair info from a given fingerprint or uid
def gpg_find_key(key_selector):
    cmd = [
        'gpg',
        '--list-keys',
        '--with-colons',
        '--fingerprint',
        '--keyid-format=long',
        key_selector
    ]
    try:
        output = subprocess.check_output(cmd, text=True)
    except FileNotFoundError as e:
        raise RuntimeError('`gpg` not found on `PATH`.') from e
    except subprocess.CalledProcessError:
        raise RuntimeError(
            'No matching key found for the provided UID/fingerprint'
        )

    key_info = {
        'fingerprint': None,
        'algorithm': None,
        'length': None,
        'curve': None,
        'uids': [],
        'chosen_uid': None
    }

    fingerprint = \
        (_normalize_fingerprint(key_selector)
         if re.fullmatch(r'[0-9A-Fa-f\s]+', key_selector or '')
         and len(_normalize_fingerprint(key_selector)) >= 32
         else None)

    # format for the output of `gpg --list-keys`
    # pub:...:<len>:<algo>:<keyid>:...
    # fpr:::::::::<PRIMARY-FINGERPRINT>::    -> fingerprint for the primary key
    # uid:::::::<Name <email>>:              -> UID(s) for the primary key
    # uid:::::::<Other Name <other@example>>:
    # sub:...:                               -> subkey (not the primary)
    # fpr:::::::::<SUBKEY-FINGERPRINT>::     -> fingerprint for the subkey
    # ...

    # this state flag tells us whether or not we're in the primary key block
    in_primary = False

    for line in output.splitlines():
        parts = line.split(':')
        tag = parts[0]
        # public key tag; the primary key info that matches
        # the given uid or fingerprint will be here
        if tag == 'pub':
            in_primary = True
            length = parts[2] or '0'
            algorithm_num = parts[3] or ''
            curve = parts[15] if len(parts) >= 16 and parts[15] else None
            key_info['length'] = int(length) if length.isdigit() else 0
            key_info['algorithm'] = \
                _PUBKEY_ALG.get(algorithm_num, f'ALG-{algorithm_num}')
            key_info['curve'] = curve
        # subkey fingerprint (if applicable)
        elif in_primary and tag == 'fpr' and key_info['fingerprint'] is None:
            normalized_fingerprint = _normalize_fingerprint(parts[9])
            if fingerprint and normalized_fingerprint != fingerprint:
                continue
            # fill in fingerprint if given keypair id was name/email
            key_info['fingerprint'] = normalized_fingerprint
        # fill in name/email from given uid
        elif in_primary and tag == 'uid':
            raw = parts[9]
            name, email = _parse_uid(raw)
            key_info['uids'].append({'raw': raw, 'name': name, 'email': email})

    if not key_info['fingerprint']:
        raise RuntimeError('Could not determine primary key fingerprint '
                           'from `gpg` output.')

    chosen_uid = None
    # identifies if Name <email> was used as the keypair identifier
    if key_selector and '<' in key_selector and '>' in key_selector:
        # since there can be multiple uids associated with a given keypair
        for uid in key_info['uids']:
            if uid['raw'] == key_selector.strip():
                chosen_uid = uid
                break
    key_info['chosen_uid'] = \
        chosen_uid or (key_info['uids'][0] if key_info['uids'] else
                       {'raw': None, 'name': None, 'email': None})

    return key_info


# helper for formatting keypair algorithm in metadata.yaml
def format_algorithm(key_info):
    algorithm = key_info.get('algorithm')
    curve = (key_info.get('curve') or '').lower()
    length = key_info.get('length') or 0
    if algorithm == 'EdDSA' and curve == 'ed25519':
        return 'Ed25519'
    elif algorithm in {'ECDSA', 'ECDH'} and key_info.get('curve'):
        return f'{algorithm}/{key_info["curve"]}'
    elif algorithm in {'RSA', 'DSA'} and length:
        return f'{algorithm}-{length}'
    else:
        return algorithm or 'unknown'
