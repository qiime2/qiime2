import os
import pkg_resources
import collections

import bibtexparser as bp


CitationRecord = collections.namedtuple('CitationRecord', ['type', 'fields'])


class Citations(dict):
    @classmethod
    def load(cls, path, package=None):
        if package is not None:
            root = pkg_resources.resource_filename(package, '.')
            root = os.path.abspath(root)
            path = os.path.join(root, path)

        with open(path) as fh:
            try:
                db = bp.load(fh)
            except Exception as e:
                raise ValueError("There was a problem loading the BiBTex file:"
                                 "%r" % path) from e
        entries = {}
        for entry in db.entries:
            id_ = entry.pop('ID')
            type_ = entry.pop('ENTRYTYPE')
            if id_ in entries:
                raise ValueError("Duplicate entry-key found in BibTex file: %r"
                                 % id_)
            entries[id_] = CitationRecord(type_, entry)

        return cls(entries)

    @classmethod
    def unformatted_citation(cls, text):
        return CitationRecord('misc', {'note': text})

    @classmethod
    def website_citation(cls, url):
        return CitationRecord('misc', {
            'note': 'No citation available. Cite plugin website.',
            'url': url
        })

    def __iter__(self):
        yield from self.values()
