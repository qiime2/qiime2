# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.core.type.util import is_visualization_type, is_collection_type


class Proxy:
    pass


class ProxyResult(Proxy):
    def __init__(self, future, selector, signature=None):
        """We have a future that represents the results of some QIIME 2 action,
        and we have a selector indicating specifically which result we want
        """
        self._future_ = future
        self._selector_ = selector
        self._signature_ = signature

    @property
    def _archiver(self):
        return self.result()._archiver

    @property
    def type(self):
        return self._archiver.type

    @property
    def uuid(self):
        return self._archiver.uuid

    @property
    def format(self):
        return self._archiver.format

    @property
    def citations(self):
        return self._archiver.citations

    def result(self):
        return self.get_element(self._future_.result())

    def get_element(self, results):
        """Get the result we want off of the future we have
        """
        return getattr(results, self._selector_)


# TODO: Can put all artifact API methods on this class because we know what
# type it will be (maybe not uuid)
class ProxyArtifact(ProxyResult):
    """This represents a future Artifact that is being returned by a Parsl app
    """
    def __repr__(self):
        if self._signature_ is None:
            return f'<artifact: Unknown Type {object.__repr__(self)}>'
        else:
            return \
                f'<artifact: {self._signature_[self._selector_].qiime_type}>'

    def view(self, type):
        """If we want to view the result we need the future to be resolved
        """
        return self.get_element(self._future_.result()).view(type)


class ProxyVisualization(ProxyResult):
    pass


class ProxyCollection(Proxy):
    def __init__(self, future, selector, signature=None):
        self._future_ = future
        self._selector_ = selector
        self._signature_ = signature

    def keys(self):
        return self.result().collection.keys()

    def values(self):
        return self.result().collection.values()

    def items(self):
        return self.result().collection.items()

    def get_element(self, results):
        """Get the result we want off of the future we have
        """
        return getattr(results, self._selector_)

    def result(self):
        return self.get_element(self._future_.result())


class ProxyResults(Proxy):
    """This represents future results that are being returned by a Parsl app
    """
    def __init__(self, future, signature):
        """We have the future results and the outputs portion of the signature
        of the action creating the results
        """
        self._future_ = future
        self._signature_ = signature

    def __iter__(self):
        """Give us a ProxyArtifact for each result in the future
        """
        for s in self._signature_:
            yield self._create_proxy(s)

    def __getattr__(self, attr):
        """Get a particular ProxyArtifact out of the future
        """
        return self._create_proxy(attr)

    def __getitem__(self, index):
        return self._create_proxy(list(self._signature_.keys())[index])

    def __repr__(self):
        lines = []
        lines.append('%s (name = value)' % self.__class__.__name__)
        lines.append('')

        max_len = -1
        for field in self._signature_:
            if len(field) > max_len:
                max_len = len(field)

        for field, value in zip(self._signature_, self):
            field_padding = ' ' * (max_len - len(field))
            lines.append('%s%s = %r' % (field, field_padding, value))

        max_len = -1
        for line in lines:
            if len(line) > max_len:
                max_len = len(line)
        lines[1] = '-' * max_len

        return '\n'.join(lines)

    def result(self):
        return self._future_.result()

    def _create_proxy(self, selector):
        qiime_type = self._signature_[selector].qiime_type

        if is_collection_type(qiime_type):
            return ProxyCollection(self._future_, selector, self._signature_)
        elif is_visualization_type(qiime_type):
            return ProxyVisualization(
                self._future_, selector, self._signature_)

        return ProxyArtifact(self._future_, selector, self._signature_)
