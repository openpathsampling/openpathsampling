__author__ = 'jan-hendrikprinz'

import openpathsampling as paths
from openpathsampling.todict import ops_object

@ops_object
class PathMoveChange(object):
    '''
    A class that described the concrete realization of a PathMove.

    Attributes
    ----------
    mover : PathMover
        The mover that generated this PathMoveChange
    generated : list of Sample
        A list of newly generated samples by this perticular move.
        Only used by node movers like RepEx or Shooters
    subchanges : list of PathMoveChanges
        the PathMoveChanges created by submovers
    details : Details
        an object that contains MoveType specific attributes and information.
        E.g. for a RandomChoiceMover which Mover was selected.
    '''

    @staticmethod
    def _indent(s):
        spl = s.split('\n')
        spl = [' |  ' + p if p[0] == ' ' else ' +- ' + p for p in spl]
        return '\n'.join(spl)

    def __init__(self, subchanges=None, samples=None, mover=None, details=None):
        self._len = None
        self._collapsed = None
        self._results = None
        self._trials = None
        self._accepted = None
        self.mover = mover
        if subchanges is None:
            self.subchanges = list()
        else:
            self.subchanges = subchanges

        if samples is None:
            self.samples = list()
        else:
            self.samples = samples
        self.details = details

    @property
    def subchange(self):
        """
        Return the single/only sub-pathmovechange if there is only one.

        Returns
        -------
        PathMoveChange
        """
        if len(self.subchanges) == 1:
            return self.subchanges[0]
        else:
            # TODO: might raise exception
            return None

    def __iter__(self):
        yield self
        for subchange in self.subchanges:
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
            return self.subchanges[item]

        if type(item) is list:
            # this is assumed to be a tree
            if item[0] is self.mover:
                if len(item) > 1:
                    for ch in self.subchanges:
                        r = ch[item[1]]
                        if r is not None:
                            return r
                    return None
                else:
                    return self
            else:
                return None


    def __reversed__(self):
        for subchange in self.subchanges:
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

    @property
    def submovers(self):
        return [ch.mover for ch in self.subchanges]

    @staticmethod
    def _check_tree(tree, branch, match):
        WILDCATS = {
            '*' : lambda s : slice(0,None),
            '.' : lambda s : slice(1,2),
            '?' : lambda s : slice(0,2),
            ':' : lambda s : slice(*map(int, s.split(':')))
        }
        MATCH_ONE = ['.', '?', '*']

#        print branch, 'in', tree

        if branch[0] not in MATCH_ONE and not match(tree[0], branch[0]):
            return False
        else:
            if len(branch) > 1:
                sub = branch[1]
                sub_branch = [branch[0]] + branch[2:]
                if type(sub) is str:
                    region = None
                    for wild in WILDCATS:
                        if wild in sub:
                            region = WILDCATS[wild](sub)
                            break

                    if region is None:
                        raise ValueError('Parse error. ONLY ' + str(WILDCATS.values()) + ' as wildcats allowed.')

                    if region.start < len(tree):
                        # check that there are enough children to match
                        for left in range(*region.indices(len(tree))):

                            sub_tree = [tree[0]] + tree[1+left:]
                            if PathMoveChange._check_tree(sub_tree, sub_branch, match):
                                return True

                    return False
                else:
                    if len(tree) > 1:
                        if not PathMoveChange._check_tree(tree[1], sub, match):
                            return False
                        else:
                            # go to next sub in branch
                            if len(branch) > 2:
                                if len(tree) > 2:
                                    return PathMoveChange._check_tree([tree[0]] + tree[2:], sub_branch, match)
                                else:
                                    return False
                    else:
                        # still branch, but no more tree
                        return False

        return True

    @staticmethod
    def _default_match(original, test):
        if isinstance(test, paths.PathMoveChange):
            return original is test
        elif isinstance(test, paths.PathMover):
            return original.mover is test
        elif issubclass(test, paths.PathMover):
            return original.mover.__class__ is test
        else:
            return False

    def _check_head_node(self, items):
        tree = self.tree()
        return self._check_tree(tree, items, self._default_match)


    def __contains__(self, item):
        """
        Check if a pathmover, pathmovechange or a tree is in self

        The tree structure is as follows

        1. A tree consists of nodes
        2. Each node can have zero, one or more children
        3. Each child is a node itself

        The tree structure in openpathsampling is expressed as

        1. The tree structure is given as a nested list.
        2. The first element in the list is the node
        3. Element 2 to N are the children.
        4. Children are always wrapped in brackets
        5. An element can be a PathMover instance or PathMoveChange instance

        node = [element, [child1], [child2], ... ]

        A tree can be a subtree if the subtree (always starting from the top)
        fits on top of the tree to match. Here child nodes are ignored as long
        as the mask of the subtree fits.

        In searching wildcats are allowed. This works

        1. slice(start, end) means an a number of arbitrary children between
            start and end-1
        2. '*' means an arbitrary number of arbitrary children. Equal to slice(0, None)
        3. None or '.' means ONE arbitrary child. Equal to slice(1,2)
        4. '?' means ONE or NONE arbitrary child. Equal to slice(0,2)

        Examples
        --------
        >>> tree1 = [mover1, [mover2]]
        >>> tree2 = [mover1, [mover2], [mover3]]
        >>> tree3 = [mover1, [mover2], [mover4, [mover5]]]
        >>> tree4 = []

        Parameters
        ----------
        item : PathMover, PathMoveChange, PathMoveTree

        """
        if isinstance(item, paths.PathMover):
            return item in self.map_post_order(lambda x : x.mover)
        elif isinstance(item, paths.PathMoveChange):
            return item in iter(self)
        elif type(item) is list:
            return self._check_head_node(item)

            # Disable checking for submoves for now

            # the head node did not fit so continue trying subnodes
#            for sub in self.subchanges:
#                if item in sub:
#                    return True

            return False

        else:
            raise ValueError('Only PathMovers or PathMoveChanges or trees ' +
                             'can be tested.')

    def __str__(self):
        return self.__class__.__name__[:-14] + '(' + str(self.idx.values()) + ')'

    def __repr__(self):
        return self.__class__.__name__[:-14] + '(' + str(self.idx.values()) + ')'

    def tree(self):
        return [self] + [ ch.tree() for ch in self.subchanges]

    def map_tree(self, fnc):
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
        return [fnc(self)] + [ ch.map_tree(fnc) for ch in self.subchanges]

    def movetree(self):
        return self.map_tree(lambda x : x.mover)

    def keytree(self, movepath=None):

        if movepath is None:
            return self.keytree([self.mover])

        result = list()
        result.append( ( movepath, self ) )
        mp = []
        for sub in self.subchanges:
            subtree = sub.keytree()
            result.extend([ ( movepath + mp + [m[0]], m[1] ) for m in subtree ])
            mp.extend([subtree[-1][0]])

        return result

    def map_post_order(self, fnc, **kwargs):
        """
        Traverse the tree of pathmovechanges in post-order applying a function

        This maps the underlying tree of pathmovechanges and applies the
        given function at each node returning a list of the results. Post-order
        will result in the order in which samples are generated. That means
        that subchanges are called first BEFORE the node itself is evaluated.

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
        that subchanges are called first BEFORE the node itself is evaluated.

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
        for mp in self.subchanges:
            output.extend(mp.level_post_order(fnc, level + 1, **kwargs))
        output.append((level, fnc(self, **kwargs)))

        return output

    def map_pre_order(self, fnc, **kwargs):
        """
        Traverse the tree of pathmovechanges in pre-order applying a function

        This maps the underlying tree of pathmovechanges and applies the
        given function at each node returning a list of the results. Pre-order
        means that subchanges are called AFTER the node itself is evaluated.

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
        means that subchanges are called AFTER the node itself is evaluated.

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

        for mp in self.subchanges:
            output.extend(mp.level_pre_order(fnc, level + 1, **kwargs))

        return output

    @property
    def collapsed_samples(self):
        """
        Return a collapsed set of samples with non used samples removed

        This is the minimum required set of samples to keep the `PathMoveChange`
        correct and allow to target sample set to be correctly created.
        These are the samples used by `.closed`

        Example
        -------
        Assume that you run 3 shooting moves for replica #1. Then only the
        last of the three really matters for the target sample_set since #1
        will be replaced by #2 which will be replaced by #3. So this function
        will return only the last sample.

        See also
        --------
        PathMoveChange.closed
        PathMoveChange.reduced()

        """
        if self._collapsed is None:
            s = paths.SampleSet([]).apply_samples(self.results)

            # keep order just for being thorough
            self._collapsed = [
                samp for samp in self.results
                if samp in s
            ]

        return self._collapsed

    @property
    def trials(self):
        """
        A list of all samples generated for this move.

        This contains accepted, rejected as well as samples that are
        redundant at the end because they we will be replaced by other
        samples already within the move

        Returns
        -------
        list of Sample
            the ordered list of all generated samples

        Notes
        -----
        E.g. Imagine a Sequential move that creates 20 trials to replace
        sample in ensemble 1. This will return all 20 samples although
        at most one sample is necessary to represent the full move, i.e. the
        last accepted one.

        """
        if self._trials is None:
            self._trials = self._get_trials()

        return self._trials

    @property
    def all_samples(self):
        """
        Returns a list of all samples generated during the PathMove.

        This includes all accepted and rejected samples (which does NOT
        include hidden samples yet)

        """
        if self._trials is None:
            self._trials = self._get_trials()
        return self._trials

    @property
    def accepted(self):
        """
        Returns if this particular move was accepted.

        Mainly used for rejected samples.
        """
        if self._accepted is None:
            self._accepted = len(self.results) > 0

        return self._accepted

    def __add__(self, other):
        return SequentialPathMoveChange([self, other])

    def apply_to(self, other):
        """
        Standard is to apply the list of samples contained
        """

        new_sample_set = paths.SampleSet(other).apply_samples(self.results)
        new_sample_set.movepath = self
        return new_sample_set

    def __call__(self, other):
        return self.apply_to(other)

    @property
    def results(self):
        """
        Returns a list of all samples that are accepted in this move

        This contains unnecessary, but accepted samples, too.
        """
        if self._results is None:
            self._results = self._get_results()

        return self._results

    @property
    def finals(self):
        """
        Returns a list of all final samples that are accepted in this move

        This contains unnecessary, but accepted samples, too.

        Notes
        -----
        E.g. in case of a Sequential move that shoots 20 times in ensemble 1.
        This will return all accepted samples which is between 0 and 20,
        although only the last accepted is relevant for the whole change.

        """
        if self._results is None:
            self._results = self._get_results()

        return self._results

    def _get_results(self):
        """
        Determines all relevant accepted samples for this move

        Includes all accepted samples also from subchanges

        Returns
        -------
        list of Sample
            the list of accepted samples for this move
        """
        return []

    def _get_trials(self):
        """
        Determines all samples for this move

        Includes all samples also from subchanges

        Returns
        -------
        list of Sample
            the list of all samples generated for this move
        """
        return []

    def __str__(self):
        if self.accepted:
            return 'SampleMove : %s : %s : %d samples' % (self.mover.cls, self.accepted, len(self.trials)) + ' ' + str(self.trials) + ''
        else:
            return 'SampleMove : %s : %s :[]' % (self.mover.cls, self.accepted)


@ops_object
class EmptyPathMoveChange(PathMoveChange):
    """
    A PathMoveChange representing no changes
    """
    def __init__(self, mover=None, details=None):
        super(EmptyPathMoveChange, self).__init__(mover=mover, details=details)

    def __str__(self):
        return ''

    @property
    def all_samples(self):
        return []

@ops_object
class SamplePathMoveChange(PathMoveChange):
    """
    A PathMoveChange representing the application of samples.

    This is the most common PathMoveChange and all other moves use this
    as leaves and on the lowest level consist only of `SamplePathMoveChange`
    """
    def __init__(self, samples, mover=None, details=None):
        super(SamplePathMoveChange, self).__init__(mover=mover, details=details)

        if type(samples) is paths.Sample:
            samples = [samples]

        self.samples = samples

    def _get_results(self):
        return [ sample for sample in self.samples if sample.accepted ]

    def _get_trials(self):
        return self.samples

@ops_object
class AcceptedSamplePathMoveChange(SamplePathMoveChange):

    def _get_trials(self):
        return self.samples

    def _get_results(self):
        return self.samples

@ops_object
class RejectedSamplePathMoveChange(SamplePathMoveChange):

    def _get_trials(self):
        return self.samples

    def _get_results(self):
        return []


@ops_object
class RandomChoicePathMoveChange(PathMoveChange):
    """
    A PathMoveChange that represents the application of a mover chosen randomly
    """
    def __init__(self, subchange, mover=None, details=None):
        super(RandomChoicePathMoveChange, self).__init__(mover=mover, details=details)
        self.subchanges = [subchange]

    def _get_results(self):
        return self.subchange.samples

    def _get_trials(self):
        return self.subchange.trials

    def apply_to(self, other):
        return self.subchange.apply_to(other)

    def __str__(self):
        return 'RandomChoice :\n' + PathMoveChange._indent(str(self.subchange))



@ops_object
class SequentialPathMoveChange(PathMoveChange):
    """
    SequentialPathMoveChange has no own samples, only inferred Sampled from the
    underlying MovePaths
    """
    def __init__(self, subchanges, mover=None, details=None):
        super(SequentialPathMoveChange, self).__init__(mover=mover, details=details)
        self.subchanges = subchanges

    def to_dict(self):
        return {
            'subchanges' : self.subchanges,
            'mover' : self.mover,
            'details' : self.details
        }

    def _get_results(self):
        samples = []
        for subchange in self.subchanges:
            samples = samples + subchange.results
        return samples

    @property
    def all_samples(self):
        samples = []
        for subchange in self.subchanges:
            samples = samples + subchange.all_samples
        return samples

    def __str__(self):
        return 'SequentialMove : %s : %d samples\n' % \
               (self.accepted, len(self.results)) + \
               PathMoveChange._indent('\n'.join(map(str, self.subchanges)))

@ops_object
class PartialAcceptanceSequentialPathMoveChange(SequentialPathMoveChange):
    """
    PartialAcceptanceSequentialMovePath has no own samples, only inferred
    Sampled from the underlying MovePaths
    """

    def _get_results(self):
        changes = []
        for subchange in self.subchanges:
            if subchange.accepted:
                changes.extend(subchange.results)
            else:
                break

        return changes

    def __str__(self):
        return 'PartialAcceptanceMove : %s : %d samples\n' % \
               (self.accepted, len(self.results)) + \
               PathMoveChange._indent('\n'.join(map(str, self.subchanges)))

@ops_object
class ConditionalSequentialPathMoveChange(SequentialPathMoveChange):
    """
    ConditionalSequentialMovePath has no own samples, only inferred Samples
    from the underlying MovePaths
    """

    def _get_results(self):
        changes = []
        for subchange in self.subchanges:
            if subchange.accepted:
                changes.extend(subchange.results)
            else:
                return []

        return changes

    def __str__(self):
        return 'ConditionalSequentialMove : %s : %d samples\n' % \
               (self.accepted, len(self.results)) + \
               PathMoveChange._indent( '\n'.join(map(str, self.subchanges)))

@ops_object
class FilterSamplesPathMoveChange(PathMoveChange):
    """
    A PathMoveChange that keeps a selection of the underlying samples
    """

    def __init__(self, subchange, mover=None, details=None):
        super(FilterSamplesPathMoveChange, self).__init__(mover=mover, details=details)
        self.subchanges = [subchange]

    def _get_results(self):
        if self.mover.use_all_samples:
            # choose all generated samples
            sample_set = self.trials
        else:
            # chose only accepted ones!
            sample_set = self.results

        # allow for negative indices to be picked, e.g. -1 is the last sample
        samples = [ idx % len(sample_set) for idx in self.mover.selected_samples]

        return samples

    def _get_trials(self):
        return self.subchange.trials

    def __str__(self):
        return 'FilterMove : pick samples [%s] from sub moves : %s : %d samples\n' % \
               (str(self.selected_samples), self.accepted, len(self.results)) + \
               PathMoveChange._indent( str(self.subchange) )

@ops_object
class KeepLastSamplePathMoveChange(PathMoveChange):
    """
    A PathMoveChange that only keeps the last generated sample.

    This is different from using `.reduced` which will only change the
    level of detail that is stored. This PathMoveChange will actually remove
    potential relevant samples and thus affect the outcome of the new
    SampleSet. To really remove samples also from storage you can use
    this PathMoveChange in combination with `.closed` or `.reduced`

    Notes
    -----
    Does the same as `FilterSamplesPathMoveChange(subchange, [-1], False)`
    """
    def __init__(self, subchange, mover=None, details=None):
        super(KeepLastSamplePathMoveChange, self).__init__(mover=mover, details=details)
        self.subchanges = [subchange]

    def _get_results(self):
        samples = self.subchange.samples
        if len(samples) > 1:
            samples = [samples[-1]]

        return samples

    def _get_trials(self):
        return self.subchange.trials

    def __str__(self):
        return 'Restrict to last sample : %s : %d samples\n' % \
               (self.accepted, len(self.results)) + \
               PathMoveChange._indent( str(self.subchange) )

@ops_object
class PathSimulatorPathMoveChange(PathMoveChange):
    """
    A PathMoveChange that just wraps a subchange and references a PathSimulator
    """

    def __init__(self, subchange, mover=None, details=None):
        super(PathSimulatorPathMoveChange, self).__init__(mover=mover, details=details)
        self.subchanges = [subchange]

    def _get_results(self):
        return self.subchange.samples

    def _get_trials(self):
        return self.subchange.trials

    def __str__(self):
        return 'PathSimulatorStep : %s : Step # %d with %d samples\n' % \
               (str(self.mover.pathsimulator.cls), self.details.step, len(self.results)) + \
               PathMoveChange._indent( str(self.subchange) )