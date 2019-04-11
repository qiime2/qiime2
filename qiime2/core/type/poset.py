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

        def __repr__(self):
            return "<%r>" % self.item

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

        def iter_descendants(self):
            seen = set()
            descendants = [self.children]
            while descendants:
                for d in descendants.pop(0):
                    if d not in seen:
                        yield d
                        seen.add(d)
                        descendants.append(d.children)

        def shared_descendants(self, other):
            descendants = set(self.iter_descendants())
            for d in other.iter_descendants():
                if d in descendants:
                    yield d

    def __init__(self, *items):
        self._maximum_antichain = []
        self._minimum_antichain = []

        for item in items:
            self.add(item)

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
                children_remaining = set(node.children) - seen
                if children_remaining:
                    new_queue.append(node)
                    continue

                if node not in seen:
                    yield node
                    seen.add(node)
                    new_queue.extend(node.parents)

            queue = new_queue

    def add(self, new):
        new = self.Node(new)

        continue_search = []
        new_max = []
        for node in self._maximum_antichain:
            greater = node <= new
            smaller = new <= node
            equal = greater and smaller

            if equal:
                return
            elif greater:
                # node not max, but wait until end so as to only add
                # new to max one time
                node.replace_parent(None, new)
            elif smaller:
                new_max.append(node)  # node is still max
                continue_search.append((node, node.children))
            else:  # neither
                new_max.append(node)  # node is still max

        if len(new_max) < len(self._maximum_antichain) or not continue_search:
            self._maximum_antichain = new_max + [new]

        # handle minimum antichain before descending arbitrarily
        new_min = []
        search_above = []
        for node in self._minimum_antichain:
            greater = node <= new
            smaller = new <= node
            equal = greater and smaller

            if equal:
                return
            elif smaller:
                new.replace_parent(None, node)
            elif greater:
                new_min.append(node)
                search_above.append((node, node.parents))
            else:
                new_min.append(node)

        if len(new_min) < len(self._minimum_antichain) or not search_above:
            self._minimum_antichain = new_min + [new]
            # still need to descend through continue_search...

        descend = [continue_search]
        visited = set()
        while descend:
            for (parent, children) in descend.pop(0):
                found = False
                for node in children:
                    if node in visited:
                        continue

                    visited.add(node)
                    greater = node <= new
                    smaller = new <= node
                    equal = greater and smaller

                    if equal:
                        return
                    elif greater:
                        found = True
                        node.replace_parent(parent, new)
                        parent.replace_child(node, new)
                    elif smaller:
                        found = True
                        descend.append([(node, node.children)])

                if not found:
                    new.replace_parent(None, parent)

        ascend = [search_above]
        while ascend:
            for (child, parents) in ascend.pop(0):
                found = False
                if child in visited:
                    continue  # child is same as node in the descend routine
                for node in parents:
                    greater = node <= new
                    smaller = new <= node
                    equal = greater and smaller

                    if equal:
                        return
                    elif greater:
                        # Don't add node to visited as we're about to look at
                        # it again
                        found = True
                        ascend.append([(node, node.parents)])
                    elif smaller:
                        visited.add(node)  # do add it to visited
                        found = True
                        child.replace_parent(node, new)
                        node.replace_child(child, new)

                if not found:
                    new.replace_child(None, child)
