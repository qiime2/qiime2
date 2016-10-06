# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from qiime.core.util import LateBindingAttribute


class Action:
    __call__ = LateBindingAttribute('_callable.__call__')
    async = LateBindingAttribute('_callable.async')

    # Private constructor
    @classmethod
    def _from_callable(cls, callable, name, description, citations, source):
        """

        Parameters
        ----------
        callable : qiime.core.callable.Callable
        name : str
            Human-readable name for this action.
        description : str
            Human-readable description for this action.
        citations : list
            List of citations
        source : str
            Markdown text defining/describing this action's computation.

        """
        self = cls.__new__(cls)

        self._callable = callable
        self.name = name
        self.description = description
        self.citations = citations
        self.source = source

        return self

    def __init__(self):
        raise NotImplementedError(
            "%s constructor is private." % self.__class__.__name__)

    @property
    def id(self):
        return self._callable.id

    @property
    def signature(self):
        return self._callable.signature

    def __repr__(self):
        return repr(self._callable)
