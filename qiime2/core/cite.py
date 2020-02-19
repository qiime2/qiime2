# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources
import collections

import bibtexparser as bp

CitationRecord = collections.namedtuple('CitationRecord', ['type', 'fields'])


class Citations(collections.OrderedDict):
    @classmethod
    def load(cls, path, package=None):
        if package is not None:
            root = pkg_resources.resource_filename(package, '.')
            root = os.path.abspath(root)
            path = os.path.join(root, path)

        parser = bp.bparser.BibTexParser()
        # Downstream tooling is much easier with unicode. For actual latex
        # users, use the modern biber backend instead of bibtex
        parser.customization = bp.customization.convert_to_unicode
        with open(path) as fh:
            try:
                db = bp.load(fh, parser=parser)
            except Exception as e:
                raise ValueError("There was a problem loading the BiBTex file:"
                                 "%r" % path) from e

        entries = collections.OrderedDict()
        for entry in db.entries:
            id_ = entry.pop('ID')
            type_ = entry.pop('ENTRYTYPE')
            if id_ in entries:
                raise ValueError("Duplicate entry-key found in BibTex file: %r"
                                 % id_)
            entries[id_] = CitationRecord(type_, entry)

        return cls(entries)

    def __iter__(self):
        return iter(self.values())

    def save(self, f):
        entries = []
        for key, citation in self.items():
            entry = citation.fields.copy()
            entry['ID'] = key
            entry['ENTRYTYPE'] = citation.type

            entries.append(entry)

        db = bp.bibdatabase.BibDatabase()
        db.entries = entries

        writer = bp.bwriter.BibTexWriter()
        writer.order_entries_by = tuple(self.keys())

        owned = False
        if type(f) is str:
            f = open(f, 'w')
            owned = True
        try:
            bp.dump(db, f, writer=writer)
        finally:
            if owned:
                f.close()
