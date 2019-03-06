
from abc import ABCMeta, abstractmethod


from qiime2.core.type.grammar import IncompleteExp, TypeExp, PredicateExp


def instantiate(cls):
    return cls()


class _BaseTemplate(metaclass=ABCMeta):
    public_proxy = set()

    @abstractmethod
    def __eq__(self, other):
        raise NotImplementedError

    def __hash__(self, other):
        return 0

    @abstractmethod
    def get_name(self):
        raise NotImplementedError

    @abstractmethod
    def get_kind(self):
        raise NotImplementedError

    @abstractmethod
    def validate_union(self, other):
        raise NotImplementedError

    @abstractmethod
    def validate_intersection(self, other):
        raise NotImplementedError

    @abstractmethod
    def is_element(self, value, expr):
        raise NotImplementedError

    def collapse_intersection(self, other):
        return NotImplemented

    def is_subtype(self, other):
        return self.get_name() == other.get_name()

    def is_supertype(self, other):
        return self.get_name() == other.get_name()


class TypeTemplate(_BaseTemplate):
    def __new__(cls, *args, **kwargs):
        self = super().__new__(cls)
        self.__init__(*args, **kwargs)

        if list(self.get_field_names()):
            return IncompleteExp(self)
        else:
            return TypeExp(self)

    def specialize(self, fields):
        return self

    @abstractmethod
    def get_field_names(self):
        raise NotImplementedError

    @abstractmethod
    def validate_fields(self, fields, expr):
        raise NotImplementedError

    @abstractmethod
    def validate_predicate(self, predicate, expr):
        raise NotImplementedError


class PredicateTemplate(_BaseTemplate):
    def __new__(cls, *args, **kwargs):
        self = super().__new__(cls)
        self.__init__(*args, **kwargs)
        return PredicateExp(self)

    @abstractmethod
    def __hash__(self, other):
        raise NotImplementedError
