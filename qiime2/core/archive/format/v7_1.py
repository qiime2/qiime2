# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.core.archive.format.v7_0 as v7_0


class ArchiveFormat(v7_0.ArchiveFormat):
    """QIIME 2 Archive Format Version 7.1

    New Features
    ------------

    Annotation[Signature]
        `Signature` has been added as a new supported Annotation sub-type.

        This sub-type acts as a mechanism for users to 'self-sign' a Result.

        A new Signature can be created and attached to a Result either via
        the Python API or the CLI.

        `annotations` directory structure (containing an example Signature):

            annotations/
            ├── uuid1/
            │   ├── metadata.yaml
            │   ├── signature.gpg
            │   ├── manifest.sha512

        The `metadata.yaml` file for Signatures will contain all standard
        Annotation fields, plus Signature-specific fields:

        Standard fields:
            `id`
            `name`
            `type`
            `created_at`
            `root_result_uuid`
            `referenced_result_uuid`

        Signature-specific fields:

            `algorithm`
                The algorithm used to create the key pair being used
                for Signature creation.

            `checksum_digest`
                The calculated sha512 checksum of the root
                checksums.sha512 file.

            `signer_name`
                Name of the user creating the Signature.

            `signer_email`
                (Optional) Email address of the user creating the Signature.

            `fingerprint`
                The calculated sha512 checksum of the public key.

        The `signature.gpg` file will contain the calculated signature
        of the private key and the root level checksums.sha512 file.
    """
    pass
