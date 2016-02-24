# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

class MetaType(type):
    def __getitem__(cls, arg):
        typetree = TypeTree(arg)
        #memoize here ideally
        return self.__class__('TypeExpr', (_Type, _TypeExpr), {'__typetree__': typetree})

class _TypeExpr(metaclass=Final):
    __typetree__ = TypeTree()

class _Type(metaclass=MetaType):
    Interfaces = tuple()

    @classmethod
    def serialize(cls, directory, subcontexts):
        pass

    @classmethod
    def deserialize(cls, directory, subcontexts):
        pass

class Type(_Type):
    pass

class ArtifactType(TypeTree):
    pass

# Some file-------------------------------------



class MySubInterface(Interface):
    def my_method(self):
        pass

class MySubInterface2(Interface):
    pass

class MyArtifactType(ArtifactType):
    Interfaces = (MySubInterface, MySubInterface2)

    @classmethod
    def serialize(cls, directory, contexts):
        int1, int2 = contexts
        int1.my_method(directory.create_file())

    @classmethod
    def deserialize(cls, directory):
        pass


class Flooble(ArtifactType):
    Interfaces = tuple(MyThingy)
