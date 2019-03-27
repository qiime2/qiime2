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
            decendents = set(self.iter_decendents())
            for d in other.iter_decendents():
                if d in decendents:
                    yield d

    def __init__(self, *items):
        self._maximum_antichain = []
        self._minimum_antichain = []

        print("{{{CREATE:::", items)
        for item in items:
            print('  {{Insert::', item)
            self.add(item)
            print('  }}After:: min: ', self.minimum_antichain, '::: max: ',  self.maximum_antichain)
        print("}}}END:::")

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
            print("    Max extended")

#            if not continue_search:
#                # check minimum
#                found = False
#                for node in self._minimum_antichain:
#                    greater = node <= new
#                    smaller = new <= node
#                    equal = greater and smaller
#
#                    assert not equal  # shouldn't happen
#                    assert not smaller
#                    if greater:
#                        # we know that the maximum antichain includes `new`
#                        # now, so there should not be anything above the
#                        # minimum which is smaller than `new`, as that would
#                        # imply a relation to one of the chains dominated by
#                        # an element max antichain which would preclude `new`
#                        # from being in the max antichain
#                        found = True
#                        new.replace_child(None, node)
#                if not found:
#                    print("    Min extended from max")
#                    self._minimum_antichain.append(new)
#            print("    Exit Simple")
#            return

        # handle minimum antichain before decending arbitrarily
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
            print("    Min extended: ", len(new_min) < len(self._minimum_antichain))
            self._minimum_antichain = new_min + [new]
            # still need to decend through continue_search...


        decend = [continue_search]
        visited = set()
        while decend:
            print("    decend:", decend)
            for (parent, children) in decend.pop(0):
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
                        print("      insert", node, ' <- ', new, ' <- ', parent)
                        node.replace_parent(parent, new)
                        parent.replace_child(node, new)
                    elif smaller:
                        found = True
                        decend.append([(node, node.children)])

                if not found:
                    print("      insert ", new, ' <- ', parent)
                    new.replace_parent(None, parent)

        ascend = [search_above]
        while ascend:
            print("    ascend:", ascend)
            for (child, parents) in ascend.pop(0):
                print("      child: ", child, " parents: ", parents)
                found = False
                if child in visited:
                    print("      skip: ", child)
                    continue  # child is the same as node in the decend routine
                for node in parents:
#                    if node in visited:
#                        print("      skipn", node)
#                        continue

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
                        print("       insert", child, ' <- ', new, ' <- ', node)
                        child.replace_parent(node, new)
                        node.replace_child(child, new)

                if not found:
                    print("      insert ", child, ' <- ', new)
                    new.replace_child(None, child)


