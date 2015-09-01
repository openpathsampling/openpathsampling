__author__ = 'jan-hendrikprinz'

import itertools
import collections
import random

class TreeSetMixin(object):
    """
    A mixin that provides basic handling for sets of trees.

    A tree is basically a node with children of the same type. The mixin
    requires to implement `.subnodes` to contain a list of children.

    The `__contains__` operator requires `_default_match` to be implemented.
    Otherwise is defaults to a comparison between two elements.

    A tree set means that it actually does not represent a single tree but a whole
    group of trees. The description works by assuming that all leaves are actually
    different choices of a single leave.

    The main difficulty is that now leaves can have two meaning and we need to
    """

    NODE_TYPE_NONE = 0
    NODE_TYPE_ALL = 1
    NODE_TYPE_ONE = 2
    NODE_TYPE_ACCUMULATE = 3
    NODE_TYPE_POWER = 4
    NODE_TYPE_CUSTOM = 5

    _node_type = NODE_TYPE_NONE

    @property
    def head(self):
        return self[0]

    @property
    def tail(self):
        return TupleTree._tail(self)

    @property
    def branches(self):
        return self._subnodes

    @property
    def _subnodes(self):
        return []

    @staticmethod
    def _tail(obj):
        if len(obj._subnodes) > 0:
            return TreeSetMixin._tail(obj[-1])

        return obj[0]

    @property
    def identifier(self):
        """
        A unique identifier to build the unique key for a position in a tree

        Returns
        -------
        hashable object
            the unique (hashable) key to identify each node

        Notes
        -----
        This is often specific to the node type and hence overridden by the
        target tree
        """
        return hex(id(self))

    @property
    def _choices(self):
        if self._node_type == self.NODE_TYPE_ALL:
            return [[sub for sub in self._subnodes]]
        elif self._node_type == self.NODE_TYPE_ACCUMULATE:
            return [self._subnodes[:n+1] for n in range(0, len(self._subnodes))]
        elif self._node_type == self.NODE_TYPE_ONE:
            leaves = [[sub] for sub in self._subnodes]

            unique_leaves = list(set([ tuple(l) for l in leaves ]))

            return unique_leaves
        elif self._node_type == self.NODE_TYPE_POWER:
            s = list(self._subnodes)
            return itertools.chain.from_iterable(
                itertools.combinations(s, r) for r in range(len(s)+1)
            )
        elif self._node_type == self.NODE_TYPE_NONE:
            return []
        elif hasattr(self._node_type, '__call__'):
            return self._node_type()
        else:
            raise RuntimeError('Node has no type !')

    def treeprint(self):
        return str(self.head) + "\n" + TreeSetMixin._indent("\n".join(map(lambda x : x.treeprint(), self.branches)))

    def locate(self, item):
        l = [key for key, value in self.locators().iteritems() if self._default_match(value, item)]
        if len(l) == 0:
            return None
        elif len(l) == 1:
            return l[0]
        else:
            return l

    def pick(self, item):
        loc = self.locate(item)
        if loc is None:
            return loc
        elif type(loc) is list:
            return [self[l] for l in loc]
        else:
            return self[loc]

    @property
    def is_sequential(self):
        return len(self._choices) < 2

    @property
    def deterministic(self):
        if not hasattr(self, '_deterministic'):
            if not self.is_sequential:
                self._deterministic = False
            else:
                self._deterministic = True
                for node in self._subnodes:
                    if not node.deterministic:
                        self._deterministic = False

        return self._deterministic

    @property
    def unique(self):
        """
        Return the smallest tree of tuples that uniquely represents this tree
        """
        ret = []
        if len(self._subnodes) == 0:
            ret = [self.identifier]
        elif self.deterministic:
            ret = [self.identifier]
        elif len(self._subnodes) == 1:
            ret = [self.identifier, self._subnodes[0].unique]
        else:
            ret = [self.identifier]
            if self.is_sequential:
                for sub in self._subnodes:
                    ret.append(sub.unique)
            else:
                for sub in self._subnodes:
                    ret.append(sub.unique)

        return TupleTree(ret)

    @property
    def enum(self):
        """
        Return a generator of all possible choices of this tree
        """
        ret = TupleTree([self.identifier])
        if self.deterministic:
            yield ret
        else:
            for leave in self._choices:
                if len(leave) == 0:
                    yield ret
                else:
                    for l in itertools.product(ret, *map(lambda x : x.enum, leave)):
                        yield TupleTree(l)

    @classmethod
    def _in_tree(cls,
                 tree,
                 test,
                 node_match_fnc,
                 leave_fnc=None,
                 leave_n = 0,
                 tree_branch_n = 0,
                 test_branch_n = 0
    ):

        if leave_fnc is None:
            leave_fnc =  lambda x : x._choices

        WILDCARDS = {
            '*' : lambda s : slice(0,None),
            '.' : lambda s : slice(1,2),
            '?' : lambda s : slice(0,2),
            ':' : lambda s : slice(*map(int, s.split(':'))),
            None: lambda s : slice(1,2)
        }
        MATCH_ONE = ['.', '?', '*']

        # print leave_n, '/', len(tree._leaves), start, '/', len(tree._leaves[leave_n]), tree.__class__.__name__, tree.identifier,  match(tree.identifier, branch[0]), branch

        if test[0] not in MATCH_ONE and not node_match_fnc(tree, test[0]):
            return False
        else:
            if len(test) + test_branch_n < 2:
                return True
            else:
                sub = test[test_branch_n+1]
                if type(sub) is str:
                    region = None
                    for wild in WILDCARDS:
                        if wild in sub:
                            region = WILDCARDS[wild](sub)
                            break

                    if region is None:
                        raise ValueError('Parse error. ONLY ' + str(WILDCARDS.values()) + ' as wildcards allowed.')

                    if leave_n < len(leave_fnc(tree)):
                        leave = leave_fnc(tree)[leave_n]
                        if region.start <= len(leave):
                            # check that there are enough children to match
                            for left in range(*region.indices(len(leave) - 1)):
                                if cls._in_tree(tree, test, node_match_fnc, leave_fnc, leave_n, tree_branch_n + left, test_branch_n + 1):
                                    return True

                else:
                    if leave_n < len(leave_fnc(tree)):
                        leave = leave_fnc(tree)[leave_n]

                        if len(leave) > tree_branch_n:
                            if cls._in_tree(leave[tree_branch_n], sub, node_match_fnc, leave_fnc):
                                # go to next sub in branch
                                if len(test) - test_branch_n < 3:
                                    return True
                                else:
                                    if len(leave) > tree_branch_n + 1:
                                        return cls._in_tree(tree, test, node_match_fnc, leave_fnc, leave_n, tree_branch_n + 1, test_branch_n + 1)

                if leave_n < len(leave_fnc(tree)) - 1:
                    if cls._in_tree(tree, test, node_match_fnc, leave_fnc, leave_n + 1):
                        return True

                return False

    def _check_head_node(self, items):
        return self._in_tree(self, items, self._default_match)

    @staticmethod
    def _default_match(original, test):
        """
        A function determining the way single nodes are matched

        This function is used to test for __contains__

        Parameters
        ----------
        original : node
            the original node to be tested
        test : node-like
            the object a node should be tested for

        Returns
        -------
        bool
            True if the original is of type test. This depends on the
            actual implementation

        Notes
        -----
        Default is to test for equality `original == test`, often we might
        also allow for testing of classes, subclasses, e.g.
        [1,2,3] as a node could be tested for list and return True
        """
        if original == test:
            return True
        else:
            return False

    def __contains__(self, item):
        """
        Check if a node or a tree is in self

        The tree structure is as follows

        1. A tree consists of nodes
        2. Each node can have zero, one or more children
        3. Each child is a node itself

        The tree structure in openpathsampling is expressed as

        1. The tree structure is given as a nested list of lists ...
        2. The first element in the list is the node
        3. Element 2 to N are the children.
        4. Children are always wrapped in brackets

        node = [element, [child1], [child2], ... ]

        A tree can be a subtree if the subtree (always starting from the top)
        fits on top of the tree to match. Here child nodes are ignored as long
        as the mask of the subtree fits.

        In searching wildcards are allowed. This works as

        1. slice(start, end) means a number of arbitrary children between
            start and end-1
        2. '*' means an arbitrary number of arbitrary children. Equal to slice(0, None)
        3. None or '.' means ONE arbitrary child. Equal to slice(1,2)
        4. '?' means ONE or NONE arbitrary child. Equal to slice(0,2)
        5. 'n:m' is equal to slice(n,m), e.g. '0:3'

        Examples
        --------
        >>> tree1 = [mover1, [mover2]]
        >>> tree2 = [mover1, [mover2], [mover3]]
        >>> tree3 = [mover1, [mover2], [mover4, [mover5]]]
        >>> tree4 = []

        Parameters
        ----------
        item : node or tree
            the node or tree to be checked

        Returns
        -------
        bool
            True if the node is in the tree or if the subtree is in the tree

        """
        if isinstance(item, tuple) or type(item) is list:
            return self._check_head_node(item)
        else:
            for x in self:
                if self._default_match(x, item):
                    return True

        return False

    def map_tree(self, fnc, **kwargs):
        """
        Apply a function to each node and return a nested tree of results

        Parameters
        ----------
        fnc : function(node, args, kwargs)
            the function run at each node node. It is given the node
            and the optional (fixed) parameters
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        tree (fnc(node, \*\*kwargs))
            nested list of the results of the map
        """
        return TupleTree([fnc(self, **kwargs)] + [ ch.map_tree(fnc, **kwargs) for ch in self._subnodes])

    def locators(self):
        """
        Return a list of key : subtree tuples

        Returns
        -------
        list of tuple(key, subtree)
            A list of all subtrees with their respective keys
        """
        path = [self.identifier]

        result = collections.OrderedDict()
        result[TupleTree(path)] = self
        excludes = []
        for leave in self._choices:
            mp = []
            for pos, sub in enumerate(leave):
                subtree = sub.locators()
                leave_id = tuple(map(lambda x : x.identifier, leave[:pos+1]))
                if leave_id not in excludes:
                    # print tuple(mp) == leave_id[:-1], tuple(mp), leave_id[:-1]
                    result.update(
                        {TupleTree(path + mp + [key]) : m for key, m in subtree.iteritems()}
                    )
                    excludes.append(leave_id)

                mp.append(TupleTree([sub.identifier]))

        return result

    def map_post_order(self, fnc, **kwargs):
        """
        Traverse the tree in post-order applying a function

        This traverses the underlying tree and applies the given function at
        each node returning a list of the results. Post-order means
        that subnodes are called BEFORE the node itself is evaluated.

        Parameters
        ----------
        fnc : function(node, kwargs)
            the function run at each node. It is given the node and the
            optional (fixed) parameters
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        list (fnc(node, \*\*kwargs))
            flattened list of the results of the map

        Notes
        -----
        This uses the same order as `reversed()`

        See also
        --------
        map_pre_order, map_post_order, level_pre_order, level_post_order
        """
        return [ fnc(node, **kwargs) for node in reversed(self) ]

    def depth_post_order(self, fnc, level=0, **kwargs):
        """
        Traverse the tree in post-order applying a function with depth

        This traverses the underlying tree and applies the given function at
        each node returning a list of the results. Post-order means
        that subnodes are called BEFORE the node itself is evaluated.

        Parameters
        ----------
        fnc : function(node, \*\*kwargs)
            the function run at each node. It is given the node and the
            optional (fixed) parameters
        level : int
            the initial level
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        list of tuple(level, func(node, \*\*kwargs))
            flattened list of tuples of results of the map. First part of
            the tuple is the level, second part is the function result.

        See also
        --------
        map_pre_order, map_post_order, level_pre_order, level_post_order
        """

        output = list()
        for mp in self._subnodes:
            output.extend(mp.depth_post_order(fnc, level + 1, **kwargs))
        output.append((level, fnc(self, **kwargs)))

        return output

    def map_pre_order(self, fnc, **kwargs):
        """
        Traverse the tree in pre-order applying a function

        This traverses the underlying tree applies the given function at
        each node returning a list of the results. Pre-order means
        that subnodes are called AFTER the node itself is evaluated.

        Parameters
        ----------
        fnc : function(node, \*\*kwargs)
            the function run at each node. It is given the node and the
            optional parameters
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        list (fnc(node, \*\*kwargs))
            flattened list of the results of the map

        Notes
        -----
        This uses the same order as `iter()`

        See also
        --------
        map_pre_order, map_post_order, level_pre_order, level_post_order
        """
        return [ fnc(node, **kwargs) for node in iter(self) ]

    def depth_pre_order(self, fnc, level=0, **kwargs):
        """
        Traverse the tree of node in pre-order applying a function

        This traverses the underlying tree applies the given function at
        each node returning a list of the results. Pre-order means
        that subnodes are called AFTER the node itself is evaluated.

        Parameters
        ----------
        fnc : function(node, \*\*kwargs)
            the function run at each node. It is given the node and the
            optional parameters
        level : int
            the initial level
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        list of tuple(level, fnc(node, \*\*kwargs))
            flattened list of tuples of results of the map. First part of
            the tuple is the level, second part is the function result.

        See also
        --------
        map_pre_order, map_post_order, level_pre_order, level_post_order
        """

        output = list()
        output.append((level, fnc(self, **kwargs)))

        for mp in self._subnodes:
            output.extend(mp.depth_pre_order(fnc, level + 1, **kwargs))

        return output

    def __iter__(self):
        """
        Traverse the whole tree in pre-order

        Returns
        -------
        Iterator
            an iterator that traverses the tree in pre-order
        """
        yield self
        for subchange in self._subnodes:
            for change in subchange:
                yield change

    def __getitem__(self, item):
        """
        Return the n-th subchange

        Returns
        -------
        PathMoveChange
            the n-th subchange if this PathMoveChange uses underlying changes
        """

        if type(item) is int:
            if item == 0:
                return self
            elif item > 0:
                return self._subnodes[item - 1]
            elif item < 0:
                return self._subnodes[item]

        if isinstance(item, tuple):
            self._last_found = None
            def match_find(original, test):
                self._last_found = original
                return self._default_match(original, test)

            find_match = TreeSetMixin._in_tree(
                self, item, match_find
            )

            if find_match:
                return self._last_found
            else:
                raise KeyError('Key %s not found in tree' % item)

        if type(item) is list:
            # this is assumed to be a tree
            if self._default_match(item[0], self):
                if len(item) > 1:
                    for ch in self._subnodes:
                        r = ch[item[1]]
                        if r is not None:
                            return r
                    return None
                else:
                    return self
            else:
                return None


    def __reversed__(self):
        """
        Traverse the whole tree in post-order

        Returns
        -------
        Iterator
            an iterator that traverses the tree in post-order
        """
        for subchange in self._subnodes:
            for change in reversed(subchange):
                yield change

        yield self

    def __len__(self):
        """
        Returns the total number of Changes mad in a move.

        Returns
        -------
        int
            the number of (Sub)PathMoveChanges in this PathMoveChange

        """
        if self._len is None:
            self._len = len(list(iter(self)))

        return self._len

    def items(self):
        return self.locators().items()

    def iteritems(self):
        return self.locators().iteritems()

    @staticmethod
    def _indent(s):
        """
        Helper function to print indented subtrees

        Parameters
        ----------
        s : str
            string representation of a tree to be indented

        Returns
        -------
        str
            the indented representation
        """
        spl = s.split('\n')
        last = [ no for no, line in enumerate(spl) if len(line) > 0 and line[0] != ' ']
        if len(last) > 0:
            last = last[-1]
        else:
            last = len(spl)
        spl = [(' |  ' + p if no <= last else '    ' + p) if p[0] == ' ' else ' +- ' + p for no, p in enumerate(spl) if len(p)]
        return '\n'.join(spl)

    def random(self):
        """
        Generate a random choice of a tree if it has multiple possibilities

        """

        this_choice = random.choice(self._choices)

        return TupleTree([self] + [ ch.random() for ch in this_choice])


class TupleTree(tuple, TreeSetMixin):

    _node_type = TreeSetMixin.NODE_TYPE_ALL

    @staticmethod
    def contains(needle, haystack):
        def cmp_fnc(x, y):
            return x[0] is y

        def leave_fnc(x):
            return [x[1:]]

        return TreeSetMixin._in_tree(haystack, needle, cmp_fnc, leave_fnc)

    def __contains__(self, item):
        return self.contains(item, self)

    def __str__(self):
        return self.treeprint()

    @property
    def _leaves(self):
        return [self[1:]]

    @property
    def _subnodes(self):
        return list(self[1:])

    @property
    def identifier(self):
        return self[0]

    def __iter__(self):
        yield self.head
        for branch in self.branches:
            for x in branch:
                yield x

    def _repr_pretty_(self, p, cycle):
            if cycle:
                p.text('(...)')
            else:
                with p.group(4, '(', ')'):
                    p.text(str(self.head))
                    if len(self._subnodes) > 0:
                        for idx, item in enumerate(self._subnodes):
                            if idx == 0 or idx < len(self._subnodes):
                                p.text(',')
                                p.breakable()
                            p.pretty(item)
                    # else:
                    #     p.text(',')
                    #     p.breakable()