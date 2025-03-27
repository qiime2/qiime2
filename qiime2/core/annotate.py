# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


class Annotation():
    """
    General base class for all annotation subtypes.
    Stuff that's the same here:
    - sender
    - receiver
    - timestamp
    - uuid
    - type
    """


class Note(Annotation):
    """
    Note subclass, inherits from Annotations
    - input: either inline text or a .txt file
    """


class Citation(Annotation):
    """
    Citation subclass, inherits from Annotations
    Should be very similar to the Citation class
    Should probably use the CitationRecord named tuple for contents
    """


class Signature(Annotation):
    """
    Signature subclass, inherits from Annotations
    Need to figure out what's going on in here
    """
