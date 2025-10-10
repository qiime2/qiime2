# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import hashlib
import os
import pathlib
import re
import subprocess
import sys
import uuid as _uuid
import yaml

from collections import OrderedDict
from datetime import datetime

# NOTE: UPDATE ME WITH EACH NEW ANNOTATION SUB-TYPE
ANNOTATION_TYPE_LIST = ['Note', 'Signature']


class Annotation():
    """General base class for all Annotation sub-classes.

    Parameters
    ----------
    name : str
        Name of the annotation.
        For each Result object, all Annotations must have a unique name.
        e.g. The same named Annotation can be attached to multiple Results,
        but each Result cannot contain multiple Annotations with the same name.

    Properties
    ----------
    uuid
        The minted uuid4 for each new Annotation that's added.
        This will be the name of each new Annotation's sub-directory within the
        `annotations` directory and is separate from any Result's uuid.

    created_at
        The minted date/time when an Annotation is created.
        Note that this is separate from when an Annotation is attached to a
        Results object, as this can occur at multiple times.

    Returns
    -------
    obj
        An instantiated Annotation of the specified sub-class.
        Note that instantiation of the Annotation base class is not supported.

    See Also
    --------
    Note

    """

    @classmethod
    def load(cls, filepath):
        """Load an Annotation.

        Parameters
        ----------
        filepath : str
            Path to load the Annotation from.

        Returns
        -------
        obj
            The instantiated Annotation sub-class.

        Raises
        ------
        ValueError
            If no `note.txt` file is found under the
            corresponding annotation directory.

        """
        # TODO: there needs to be a better way to deal with these shared attrs
        with open(os.path.join(filepath, 'metadata.yaml'), 'r') as fh:
            meta_yaml = yaml.safe_load(fh)
            annotation_type = meta_yaml['type']

            # NOTE
            if annotation_type == 'Note':
                annotation = Note.__new__(Note)
                # Now attach Note attrs from metadata.yaml
                annotation.id = meta_yaml['id']
                annotation.name = meta_yaml['name']
                annotation.annotation_type = meta_yaml['type']
                annotation.created_at = meta_yaml['created_at']

                # Validate that `note.txt` exists
                note_fp = os.path.join(filepath, 'note.txt')
                if not os.path.exists(note_fp):
                    raise ValueError(
                        'Unable to load malformed Note with name: '
                        f'"{annotation.name}" due to missing `note.txt` file.'
                    )
                # Attach contents to Note
                else:
                    with open(note_fp, 'r') as fh:
                        annotation.contents = fh.read()

            # SIGNATURE
            elif annotation_type == 'Signature':
                annotation = Signature.__new__(Signature)
                # Now attach attrs from metadata.yaml
                annotation.id = meta_yaml['id']
                annotation.name = meta_yaml['name']
                annotation.annotation_type = meta_yaml['type']
                annotation.created_at = meta_yaml['created_at']
                # signature-specific attrs
                annotation.algorithm = meta_yaml['algorithm']
                annotation.checksum_digest = meta_yaml['checksum_digest']
                annotation.signer_name = meta_yaml['signer_name']
                annotation.signer_email = meta_yaml['signer_email']
                annotation.fingerprint = meta_yaml.get('fingerprint')

                # Validate that `signature.gpg` exists
                sig_fp = os.path.join(filepath, 'signature.gpg')
                if not os.path.exists(sig_fp):
                    raise ValueError(
                        'Unable to load malformed Signature with name: '
                        f'"{annotation.name}" due to missing '
                        '`signature.gpg` file.'
                    )

            else:
                annotation = UnknownAnnotation.__new__(UnknownAnnotation)
                annotation.id = meta_yaml['id']
                annotation.name = meta_yaml['name']
                annotation.annotation_type = meta_yaml['type']
                annotation.created_at = meta_yaml['created_at']

        return annotation

    # We never expect this to be hit as the base class for Annotations
    # shouldn't ever be instantiated - only the supported sub-classes.
    @property
    def annotation_type(cls):
        raise NotImplementedError

    def __init__(self, name):
        """
        Construction for an initialized Annotation.

        Attributes
        ----------
        name
            The user-provided name of the Annotation.

        id
            The uuid4 ID associated with the Annotation.

        created_at
            The datetime when the Annotation was created.

        Raises
        ------
        TypeError
            If the Annotation base class is instantiated.

        """
        if type(self) is Annotation:
            raise TypeError('Annotation is an abstract class'
                            ' and cannot be instantiated directly.')
        self.name = name
        self.id = _uuid.uuid4()
        self.created_at = datetime.now()

    def validate_name(self, name):
        """Validates that the given name is a valid Python idenitifier with the
        exception that `-` is allowed.

        Parameters
        ----------
        name : str
            The name to validate.

        Raises
        ------
        ValueError
            If the name passed in is not a valid Python identifier.
        """
        validate_name = name.replace('-', '_')
        if not validate_name.isidentifier():
            raise ValueError(f'Name "{name}" is not a valid Python identifier.'
                             ' Keys may contain `-` characters but must'
                             ' otherwise be valid Python identifiers. Python'
                             ' identifier rules may be found here'
                             ' https://www.askpython.com/python/'
                             'python-identifiers-rules-best-practices')

    def _write_meta_yaml(self, annotations_dir,
                         root_result_uuid, referenced_result_uuid,
                         algorithm=None, checksum_digest=None,
                         signer_name=None, signer_email=None,
                         fingerprint=None):
        """Write the contents of `metadata.yaml` for a given Annotation.

        Parameters
        ----------
        annotations_dir
            The path to the `annotations` directory within a Result object.
            Located under `provenance`.

        root_result_uuid
            The uuid of the Result object where an Annotation is being added.

        referenced_result_uuid
            The uuid of the Result object that an Annotation is referring to.
            Note that in 7.0, `root_result_uuid` and `referenced_result_uuid`
            are the same (i.e. Annotations can only refer to the same Result
            they are being attached to) but separate root and referenced uuids
            will be supported in future versions.

        Returns
        -------
        str
            The filepath where the Annotation's uuid-specific subdirectory
            containing the `metadata.yaml` file was written to.

        """
        # create the annotation directory if it's not already present
        # we don't care if it contains anything, just whether or not it exists
        if not os.path.exists(annotations_dir):
            os.mkdir(annotations_dir)

        # create the unique dir for a particular annotation
        annotation_uuid_dirname = \
            os.path.join(annotations_dir, str(self.id))
        os.mkdir(annotation_uuid_dirname)

        metadata = OrderedDict()
        metadata['id'] = str(self.id)
        metadata['name'] = self.name
        metadata['type'] = self.annotation_type
        metadata['created_at'] = self.created_at
        metadata['root_result_uuid'] = root_result_uuid
        metadata['referenced_result_uuid'] = referenced_result_uuid

        if self.annotation_type == 'Signature':
            metadata['algorithm'] = algorithm
            metadata['checksum_digest'] = checksum_digest
            metadata['signer_name'] = signer_name
            metadata['signer_email'] = signer_email
            metadata['fingerprint'] = fingerprint

        meta_yaml = os.path.join(annotation_uuid_dirname, 'metadata.yaml')
        with open(meta_yaml, 'w') as fh:
            fh.write(yaml.dump(metadata))

        return annotation_uuid_dirname


class UnknownAnnotation(Annotation):
    """Utility sub-class that handles loading newer Annotation types on an
    older version of QIIME 2 that supports Annotations.
    """
    def __init__(*args):
        raise NotImplementedError('`UnknownAnnotation` is an abstract class'
                                  ' used for handling Annotations associated'
                                  ' with future versions of QIIME 2.'
                                  ' It should not be instantiated directly.')

    _write_meta_yaml = __init__
    _write = __init__


class Note(Annotation):
    """Note sub-class, inherits from Annotations.

    Parameters
    ----------
    text : str
        Inline text that will be written inside the Note's `note.txt` file.
        This parameter is optional, but either `text` OR `filepath` must be
        provided.

    filepath : str
        Path to a file whose contents should be written inside the Note's
        `note.txt` file.
        This parameter is optional, but either `text` OR `filepath` must be
        provided.

    Properties
    ----------
    type : Note
        The type of Annotation being instantiated.

    Returns
    -------
    Note : obj
        The instantiated Note.

    See Also
    --------
    Annotation

    """
    annotation_type = 'Note'

    # NOTE: in future versions, name will become optional & the default value
    # will be the annotation's UUID (if name isn't provided by the user)
    def __init__(self, name, *, text=None, filepath=None):
        self.validate_name(name)
        # Ensure exactly one of text or filepath is provided
        if text and filepath:
            raise ValueError(
                'Cannot set both `text` and `filepath` params. '
                'Please provide either inline text or a filepath only.'
            )
        if not text and not filepath:
            raise ValueError(
                'No inputs provided to either `text` or `filepath`. '
                'Please provide either inline text or a filepath.'
            )

        if text is not None:
            self.contents = text
            self._filepath = None
        else:
            if not isinstance(filepath, (str, os.PathLike)):
                raise TypeError(
                    f'Unexpected input for `filepath`: {filepath} '
                    '`filepath` should either be a '
                    '`string` or an `os.PathLike` object.'
                )
            if not os.path.exists(filepath):
                raise ValueError(
                    f'File not found from provided `filepath`: {filepath} '
                    'Double check that the provided file exists '
                    ' and is in the expected location.'
                )

            self._filepath = str(filepath)
            self.contents = None

        # Construct Annotation class
        super().__init__(name)

    def _write(self, annotations_dir, root_result_uuid,
               referenced_result_uuid):
        """Write the contents of an instantiated Note.

        Parameters
        ----------
        annotations_dir
            The path to the `annotations` directory within a Result object.
            Located under `provenance`.

        root_result_uuid
            The uuid of the Result object where an Annotation is being added.

        referenced_result_uuid
            The uuid of the Result object that an Annotation is referring to.
            Note that in 7.0, `root_result_uuid` and `referenced_result_uuid`
            are the same (i.e. Annotations can only refer to the same Result
            they are being attached to) but separate root and referenced uuids
            will be supported in future versions.

        See Also
        --------
        _write_meta_yaml

        """
        # call _write_meta_yaml to write the stuff that's the same
        # across both input types
        annotation_uuid_dirname = \
            self._write_meta_yaml(annotations_dir,
                                  root_result_uuid,
                                  referenced_result_uuid)

        note_path = os.path.join(annotation_uuid_dirname, 'note.txt')

        if self._filepath:
            with open(self._filepath, 'rb') as fh:
                contents = fh.read()
        else:
            contents = self.contents.encode('utf-8')

        # validation for max size and parsability
        max_size = 10 * 1024 * 1024

        if len(contents) > max_size:
            raise ValueError('Note contents exceed maximum size of 10 MiB')

        try:
            contents.decode('utf-8')
        except UnicodeDecodeError:
            raise ValueError('Note contents are not valid UTF-8')

        with open(note_path, 'wb') as fh:
            fh.write(contents)


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


# SIGNATURE HELPERS
# helper for locating root_fp for a given Result
def _find_root_fp(annotations_dir, root_result_uuid):
    split_fp = annotations_dir.split(os.sep)
    root_result_uuid_index = split_fp.index(root_result_uuid)
    root_fp = os.sep.join(split_fp[0: root_result_uuid_index + 1])
    return root_fp


# helper for calculating the root level checksum digest
def _sha512_file_hex(path):
    hex = hashlib.sha512()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1024*1024), b""):
            hex.update(chunk)
    return hex.hexdigest()


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


# helper for pulling keypair info from a given fingerprint or uid
def _gpg_find_key(key_selector):
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
def _format_algorithm(key_info):
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


class Signature(Annotation):
    """Signature sub-class, inherits from Annotations.

    Creates a cryptographic signature over the Result's root checksums.sha512
    file using credentials for an existing key pair via GnuPG.

    Parameters
    ----------
    name : str
        Annotation name (validated like other Annotations).

    signer_name : str, optional
        Name associated with the key pair in GnuPG
        that will be used for signing.
        At least one of name/email must be provided for Signature creation.

    signer_email : str, optional
        Email associated with the key pair in GnuPG
        that will be used for signing.
        At least one of name/email must be provided for Signature creation.

    Returns
    -------
    Signature : obj
        The instantiated Signature.

    See Also
    --------
    Annotation
    """
    annotation_type = 'Signature'

    # NOTE: in future versions, name will become optional & the default value
    # will be the annotation's UUID (if name isn't provided by the user)
    def __init__(self, name, *, signer_uid=None, fingerprint=None):
        self.validate_name(name)

        # Ensure at least one of signer uid (name/email address) is provided
        if not signer_uid and not fingerprint:
            raise ValueError(
                'No inputs provided to either `signer_uid` or `fingerprint`. '
                'Please provide either signer ID (name and email) or '
                'fingerprint for key pair identification.'
            )

        # Construct Annotation class
        self.signer_uid = signer_uid
        self.fingerprint = fingerprint
        super().__init__(name)

    def _write(self, annotations_dir, root_result_uuid,
               referenced_result_uuid):
        """Write the contents of an instantiated Signature.

        Parameters
        ----------
        annotations_dir
            The path to the `annotations` directory within a Result object.
            Located under `provenance`.

        root_result_uuid
            The uuid of the Result object where an Annotation is being added.

        referenced_result_uuid
            The uuid of the Result object that an Annotation is referring to.
            Note that in 7.0 & 7.1, `root_result_uuid` and
            `referenced_result_uuid` are the same (i.e. Annotations can only
            refer to the same Result they are being attached to) but separate
            root and referenced uuids will be supported in future versions.
        """
        root_fp = _find_root_fp(annotations_dir, root_result_uuid)

        checksums_fp = root_fp / 'checksums.sha512'
        if not checksums_fp.exists():
            raise ValueError(
                'Unable to create Signature due to malformed Result: '
                f'missing root checksums file at "{checksums_fp}".'
            )
        checksum_digest = _sha512_file_hex(checksums_fp)

        keypair_id = self.signer_uid or self.fingerprint
        if not keypair_id:
            raise ValueError(
                'No signer identity available. `signer_uid` '
                '(e.g. "Name <email>") or `fingerprint` must be available '
                'within gpg when constructing Signature.'
            )

        key_info = _gpg_find_key(keypair_id)
        algorithm = _format_algorithm(key_info)
        fingerprint = key_info['fingerprint']
        signer_name = key_info['chosen_uid']['name']
        signer_email = key_info['chosen_uid']['email']

        annotation_uuid_dirname = \
            self._write_meta_yaml(annotations_dir, root_result_uuid,
                                  referenced_result_uuid, algorithm,
                                  checksum_digest, signer_name, signer_email,
                                  fingerprint)

        signature_dir = pathlib.Path(annotation_uuid_dirname)
        sig_fp = signature_dir / 'signature.gpg'

        env = os.environ.copy()
        # Apparently this is helpful on Unix for GPG to find
        # the correct terminal
        try:
            if sys.stdin and sys.stdin.isatty():
                env.setdefault('GPG_TTY', os.ttyname(sys.stdin.fileno()))
        except Exception:
            pass

        cmd = [
            'gpg',
            '--local-user', str(keypair_id),
            '--output', str(sig_fp),
            '--detach-sign', str(checksums_fp)
        ]

        try:
            subprocess.run(cmd, check=True, env=env)
        except FileNotFoundError as e:
            raise RuntimeError(
                'GnuPG (`gpg`) is not installed or not on `PATH`. '
                'Install GnuPG and ensure your key pair is available.'
            ) from e
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                '`gpg` signing failed. Ensure that the selected key exists '
                'and is unlocked (or that pinentry can prompt for the '
                'private key password).'
            ) from e

        if not sig_fp.exists() or sig_fp.stat().st_size == 0:
            raise RuntimeError(
                '`gpg` reported success but no signature was written at '
                f'{sig_fp!s}.'
            )

        # smoke check to ensure signature writing was successful
        try:
            subprocess.run(
                ['gpg', '--verify', str(sig_fp), str(checksums_fp)],
                check=True, env=env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
        except subprocess.CalledProcessError:
            raise RuntimeError(
                'Wrote signature.gpg but `gpg --verify` failed; '
                'signature may be invalid.'
            )
