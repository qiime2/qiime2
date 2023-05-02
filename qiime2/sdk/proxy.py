# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.core.type.util import is_visualization_type, is_collection_type
from qiime2.core.type.collection import Collection


class Proxy:
    """Base class to indicate that a given class that inherits from it is a
    proxy. Also implements some generic functionality
    """
    def __eq__(self, other):
        return self.result == other.result()

    def __ne__(self, other):
        return not (self == other)

    def _create_proxy(self, selector):
        qiime_type = self._signature_[selector].qiime_type

        if is_collection_type(qiime_type):
            return ProxyResultCollection(
                self._future_, selector, self._signature_)
        elif is_visualization_type(qiime_type):
            return ProxyVisualization(
                self._future_, selector, self._signature_)

        return ProxyArtifact(self._future_, selector, self._signature_)


class ProxyResult(Proxy):
    def __init__(self, future, selector, signature=None):
        """We have a future that represents the results of some QIIME 2 action,
        and we have a selector indicating specifically which result we want
        """
        self._future_ = future
        self._selector_ = selector
        self._signature_ = signature

    def __repr__(self):
        if self._signature_ is None:
            return f'<{self.__class__.__name__.__lower__}: Unknown Type ' \
                   f'{object.__repr__(self)}>'
        else:
            return f'<{self.__class__.__name__.__lower__}: {self.type}>'

    @property
    def _archiver(self):
        return self.result()._archiver

    @property
    def type(self):
        if self._signature_ is not None:
            return self._signature_[self._selector_].qiime_type

        return self.result.type

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
        return self._get_element_(self._future_.result())

    def _get_element_(self, results):
        """Get the result we want off of the future we have
        """
        return getattr(results, self._selector_)

    def save(self, filepath, ext=None):
        """Blocks then calls save on the result.
        """
        return self.result().save(filepath, ext=ext)


class ProxyArtifact(ProxyResult):
    """This represents a future Artifact that is being returned by a Parsl app
    """
    def view(self, type):
        """If we want to view the result we need the future to be resolved
        """
        return self._get_element_(self._future_.result()).view(type)


class ProxyVisualization(ProxyResult):
    """This represents a future Visualization that is being returned by a Parsl
       app
    """
    pass


class ProxyResultCollection(Proxy):
    def __init__(self, future, selector, signature=None):
        self._future_ = future
        self._selector_ = selector
        self._signature_ = signature

    def __len__(self):
        return len(self.collection)

    def __iter__(self):
        yield self.collection.__iter__()

    def __setitem__(self, key, item):
        self.collection[key] = item

    def __getitem__(self, key):
        return self.collection[key]

    def __repr__(self):
        return f"<{self.__class__.__name__.lower()}: {self.type}>"

    @property
    def type(self):
        # I'm note a huge fan of the fact that this may or may not need to
        # to block. If this is a return from an action (which it basically
        # always will be) we don't need to block for type. Otherwise we do.
        if self._signature_ is not None:
            return Collection[self._signature_[self._selector_]]

        return self.result().type

    @property
    def collection(self):
        return self.result().collection

    def keys(self):
        return self.collection.keys()

    def values(self):
        return self.collection.values()

    def items(self):
        return self.collection.items()

    def _get_element_(self, results):
        """Get the result we want off of the future we have
        """
        return getattr(results, self._selector_)

    def result(self):
        return self._get_element_(self._future_.result())

    def save(self, directory):
        """Blocks then calls save on the result.
        """
        return self.result().save(directory)


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

    def __eq__(self, other):
        """ Overriding the one on Proxy because we have _result not result
        """
        return self._result() == other._result()

    def _result(self):
        """ If you are calling an action in a try-except block in a pipeline,
            you need to call this method on the Results object returned by the
            action.

            This is because if the Pipeline was executed with parsl, we need to
            block on the action in the try-except to ensure we get the result
            and raise the potential exception while we are still inside of the
            try-except. Otherwise we would just get the exception whenever the
            future resolved which would likely be outside of the try-except, so
            the exception would be raised and not caught.

            If you call an action in the Python API using parsl inside of a
            context manager (a withed in Cache for instance) you also must call
            this method there to ensure you get you don't start using a
            different cache/pool/whatever before your future resolves.
        """
        return self._future_.result()
