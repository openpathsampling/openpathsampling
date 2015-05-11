__author__ = 'jan-hendrikprinz'

class TreeMixin(object):
    """
    A mixin that provides basic tree handling.

    A tree is basically a node with children of the same type. The mixin
    requires to implement
    """
    
    @property
    def subnodes(self):
        return []

    @property
    def subnode(self):
        """
        Return the single/only sub-pathmovechange if there is only one.

        Returns
        -------
        PathMoveChange
        """
        if len(self.subnodes) == 1:
            return self.subnodes[0]
        else:
            # TODO: might raise exception
            return None

    def __iter__(self):
        yield self
        for subchange in self.subnodes:
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
            return self.subnodes[item]

    def __reversed__(self):
        for subchange in self.subnodes:
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

    def key(self, change):
        tree = self.keytree()
        return [leave for leave in tree if leave[1] is change ][0][0]

    def _check_head_node(self, items):
        if isinstance(items[0], paths.PathMover):
            # a subtree of pathmovers
            if self.mover is items[0]:
                #print 'found head'
                # found current head node, check, if children match in order
                left = 0
                submovers = [ch.mover for ch in self.subnodes]
                subvalues = items[1]
                if type(subvalues) is not list:
                    subvalues = [subvalues]

                for sub in zip(subvalues[0::2], subvalues[1::2]):
                    if left >= len(self.subnodes):
                        # no more subnodes to match
                        return False
                    if sub is None:
                        # None is a placeholder so move token +1
                        left = left + 1
                    if type(sub) is dict:
                        if sub[0] is None:
                            while left < len(self.subnodes):
                                if not [self.subnodes[left].mover, [sub[1]]] in self.subnodes[left]:
                                    left = left + 1
                                else:
                                    left = left + 1
                                    break

                            if left == len(self.subnodes):
                                return False
                        elif sub[0] not in submovers[left:]:
                            #print 'missing sub', sub.keys()[0], 'in', submovers[left:]
                            return False
                        else:
                            idx = submovers.index(sub[0])
                            left = idx + 1
                            if not [sub[0], [sub[1]]] in self.subnodes[idx]:
                                #print 'try', {sub.keys()[0] : sub.values()[0]}
                                return False

                    elif isinstance(sub, paths.PathMover):
                        if sub not in submovers[left:]:
                            return False
                        idx = submovers.index(sub)
                        left = idx + 1

                return True

        elif items[0] is None or len(items) == 0:
            # means empty tree and since nothing is in every tree return true
            return True

    def __contains__(self, item):
        """
        Check if a pathmover, pathmovechange or a tree is in self

        A node is either None or a PathMover

        1. subnodes are given using a dict { parent : child }
        2. Several subnodes are given in a list. [child1, child2]
        3. A single subchange can be given as a list of length 1 or a single mover.
        4. None is a wildcat and matches everything

        Examples
        --------
        >>> tree1 = {mover1 : mover2}
        >>> tree2 = {mover1 : [mover2, mover3]}
        >>> tree3 = {mover1 : [mover2, {mover4 : [mover5]}] }
        >>> tree4 = {}

        Notes
        -----
        TODO: Add other types of nodes. e.g. explicit PathMoveChange,
        Boolean for .accepted

        Parameters
        ----------
        item : PathMover, PathMoveChange, PathMoveTree

        """
        if isinstance(item, paths.PathMover):
            return item in self.map_post_order(lambda x : x.mover)
        elif isinstance(item, paths.PathMoveChange):
            return item in iter(self)
        elif type(item) is list:
            if self._check_head_node(item):
                return True

            # Disable checking for submoves for now

            # the head node did not fit so continue trying subnodes
#            for sub in self.subnodes:
#                if item in sub:
#                    return True

            return False

        else:
            raise ValueError('Only PathMovers or PathMoveChanges can be tested.')

    def tree(self):
        return {self : [ ch.tree() for ch in self.subnodes] }

    def movetree(self):
        return {self.mover : [ ch.movetree() for ch in self.subnodes] }

    def keytree(self, movepath=None):

        if movepath is None:
            movepath = [self.mover]

        result = list()
        result.append( ( movepath, self ) )
        mp = []
        for sub in self.subnodes:
            subtree = sub.keytree()
            result.extend([ ( movepath + [mp + m[0]], m[1] ) for m in subtree ])
#            print subtree[-1][0]
            mp = mp + subtree[-1][0]


        return result

    def map_tree(self, fnc, **kwargs):
        """
        Apply a function to each node and return the tree

        Parameters
        ----------
        fnc : function(pathmovechange, args, kwargs)
            the function run at each pathmovechange node. It is given the node
            and the optional (fixed) parameters
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        tree (fnc(node, **kwargs))
            nested list of the results of the map
        """

        if len(self.subnodes) > 1:
            return { fnc(self, **kwargs) : [node.map_tree(fnc, **kwargs) for node in self.subnodes]}
        elif len(self.subnodes) == 1:
            return { fnc(self, **kwargs) : self.subnodes[0].map_tree(fnc, **kwargs)}
        else:
            return fnc(self, **kwargs)

    def map_post_order(self, fnc, **kwargs):
        """
        Traverse the tree of pathmovechanges in post-order applying a function

        This maps the underlying tree of pathmovechanges and applies the
        given function at each node returning a list of the results. Post-order
        will result in the order in which samples are generated. That means
        that subnodes are called first BEFORE the node itself is evaluated.

        Parameters
        ----------
        fnc : function(pathmovechange, args, kwargs)
            the function run at each pathmovechange node. It is given the node
            and the optional (fixed) parameters
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        list (fnc(node, **kwargs))
            flattened list of the results of the map

        Notes
        -----
        This uses the same order as `reversed()`

        See also
        --------
        map_pre_order, map_post_order, level_pre_order, level_post_order
        """
        return [ fnc(node, **kwargs) for node in reversed(self) ]

    def level_post_order(self, fnc, level=0, **kwargs):
        """
        Traverse the tree of pathmovechanges in post-order applying a function

        This maps the underlying tree of pathmovechanges and applies the
        given function at each node returning a list of the results. Post-order
        will result in the order in which samples are generated. That means
        that subnodes are called first BEFORE the node itself is evaluated.

        Parameters
        ----------
        fnc : function(pathmovechange, args, kwargs)
            the function run at each pathmovechange node. It is given the node
            and the optional parameters
        level : int
            the initial level
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        list of tuple(level, func(node, **kwargs))
            flattened list of tuples of results of the map. First part of
            the tuple is the level, second part is the function result.

        See also
        --------
        map_pre_order, map_post_order, level_pre_order, level_post_order
        """

        output = list()
        for mp in self.subnodes:
            output.extend(mp.level_post_order(fnc, level + 1, **kwargs))
        output.append((level, fnc(self, **kwargs)))

        return output

    def map_pre_order(self, fnc, **kwargs):
        """
        Traverse the tree of pathmovechanges in pre-order applying a function

        This maps the underlying tree of pathmovechanges and applies the
        given function at each node returning a list of the results. Pre-order
        means that subnodes are called AFTER the node itself is evaluated.

        Parameters
        ----------
        fnc : function(pathmovechange, args, kwargs)
            the function run at each pathmovechange node. It is given the node
            and the optional parameters
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        list (fnc(node, **kwargs))
            flattened list of the results of the map

        Notes
        -----
        This uses the same order as `iter()`

        See also
        --------
        map_pre_order, map_post_order, level_pre_order, level_post_order
        """
        return [ fnc(node, **kwargs) for node in iter(self) ]

    def level_pre_order(self, fnc, level=0, **kwargs):
        """
        Traverse the tree of pathmovechanges in pre-order applying a function

        This maps the underlying tree of pathmovechanges and applies the
        given function at each node returning a list of the results. Pre-order
        means that subnodes are called AFTER the node itself is evaluated.

        Parameters
        ----------
        fnc : function(pathmovechange, args, kwargs)
            the function run at each pathmovechange node. It is given the node
            and the optional parameters
        level : int
            the initial level
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        list of tuple(level, fnc(node, **kwargs))
            flattened list of tuples of results of the map. First part of
            the tuple is the level, second part is the function result.


        See also
        --------
        map_pre_order, map_post_order, level_pre_order, level_post_order
        """

        output = list()
        output.append((level, fnc(self, **kwargs)))

        for mp in self.subnodes:
            output.extend(mp.level_pre_order(fnc, level + 1, **kwargs))

        return output