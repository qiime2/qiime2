# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def maximum_antichain(*types):
    return POSet(*types).maximum_antichain


def minimum_antichain(*types):
    return POSet(*types).minimum_antichain


class POSet:
    class Node:
        __le_cache = {}

        def __init__(self, item):
            self.item = item
            self._children = {}
            self._parents = {}

        @property
        def children(self):
            return list(self._children)

        @property
        def parents(self):
            return list(self._parents)

        def __le__(self, other):
            a, b = comp = self.item, other.item
            if comp in self.__le_cache:
                return self.__le_cache[comp]

            r = a <= b
            self.__le_cache[comp] = r
            return r

        def __lt__(self, other):
            return self <= other and not other <= self

        def replace_parent(self, old, new):
            new._children[self] = None
            if old is not None:
                del self._parents[old]
            self._parents[new] = None

        def replace_child(self, old, new):
            new._parents[self] = None
            if old is not None:
                del self._children[old]
            self._children[new] = None

        def iter_decendents(self):
            seen = set()
            decendents = [self.children]
            while decendents:
                for d in decendents.pop(0):
                    if d not in seen:
                        yield d
                        seen.add(d)
                        decendents.append(d.children)

        def shared_decendents(self, other):
            #print([x.item for x in self.iter_decendents()])
            #print([x.item for x in other.iter_decendents()])
            decendents = set(self.iter_decendents())
            for d in other.iter_decendents():
                if d in decendents:
                    yield d

    def __init__(self, *items):
        self._maximum_antichain = []
        self._minimum_antichain = []

        print(":::START:::", items)
        for item in items:
            self.add(item)
            print(item, ' min: ', self.minimum_antichain, '::: max: ',  self.maximum_antichain)
        print(":::END:::")

    @property
    def maximum_antichain(self):
        return [n.item for n in self._maximum_antichain]

    @property
    def minimum_antichain(self):
        return [n.item for n in self._minimum_antichain]

    def __iter__(self):
        for node in self.iter_nodes():
            yield node.item

    def iter_nodes(self):
        seen = set()
        queue = self._minimum_antichain

        while queue:
            new_queue = []
            for node in queue:
                if node not in seen:
                    yield node
                    seen.add(node)
                    new_queue.extend(node.parents)
            queue = new_queue

    def add(self, new):
        new = self.Node(new)

        new_maximum_antichain = []
        smaller_than = []
        for node in self._maximum_antichain:
            if node <= new:
                if new <= node:
                    return
                else:  # node < new
                    node.replace_parent(None, new)
            else:  # node > new or node !~ new
                new_maximum_antichain.append(node)
                if new <= node:
                    smaller_than.append(node)

        if (len(new_maximum_antichain) < len(self._maximum_antichain)
                or not smaller_than):
            self._maximum_antichain = new_maximum_antichain + [new]
            if not new.children:
                # TODO: THIS IS WRONG!!!!!!!!
                # JUST BECAUSE WE ARE DISJOINT WRT MAX DOESNT MEAN WE
                # DONT HAVE A RELATION TO THE MIN
                self._minimum_antichain.append(new)
            return

        # The new entry must be below the maximum antichain.
        # Ruling out the minimum antichain makes maintaining the minimum
        # antichain a little bit easier, so we'll hold onto "smaller_than"
        # for later

        new_minimum_antichain = []
        for node in self._minimum_antichain:
            if new <= node:
                if node <= new:
                    return
                else:  # new < node
                    new.replace_parent(None, node)
            elif node <= new:
                new_minimum_antichain.append(node)

        if len(new_minimum_antichain) < len(self._minimum_antichain):
            # if there was no relation to the minimum antichain, then there
            # would be no relation to the maximum antichain, which means
            # this new entry would have been detected above during the
            # maximum antichain phase, and the minimum would have been updated
            self._minimum_antichain = new_minimum_antichain + [new]
            return

        # The new entry is not in the minimum or maximum antichains.
        # No further maintenance is needed on those items, now we just need
        # to insert it into the graph wherever it is relevant
        decend = [smaller_than]
        while decend:
            for parent in decend.pop(0):
                inserted = False
                for child in list(parent.children):
                    if child <= new:
                        if new <= child:
                            return
                        else:  # child < new and new < node
                            inserted = True
                            parent.replace_child(child, new)
                            child.replace_parent(parent, new)
                    elif new <= child:  # new < child
                        decend.append(child.children)
                if not inserted:
                    new.replace_parent(None, parent)
