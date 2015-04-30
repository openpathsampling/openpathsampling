from openpathsampling import ops_object

__author__ = 'jan-hendrikprinz'

import openpathsampling as paths
from openpathsampling.todict import ops_object

@ops_object
class PathMoveChange(object):
    #TODO: this docstring is false
    '''
    PathMoveChange is essentially a list of samples, with a few conveniences.  It
    can be treated as a list of samples (using, e.g., .append), or as a
    dictionary of ensembles mapping to a list of samples, or as a dictionary
    of replica IDs to samples. Any type is allowed as a replica ID except
    Sample or Ensemble.

    The dictionaries ensemble_dict and replica_dict are conveniences which
    should be kept consistent by any method which modifies the container.
    They do not need to be stored.

    Notes
    -----
    Current implementation is as an unordered set. Therefore we don't
    have some of the convenient tools in Python sequences (e.g.,
    slices). On the other hand, I'm not sure whether that is meaningful
    here.

    Attributes
    ----------
    samples : list of Sample
        The samples included in this set.
    ensemble_dict : dict
        A dictionary with Ensemble objects as keys and lists of Samples as
        values.
    replica_dict : dict
        A dictionary with replica IDs as keys and lists of Samples as values
    '''

    @staticmethod
    def _indent(s):
        spl = s.split('\n')
        spl = [' |  ' + p if p[0] == ' ' else ' +- ' + p for p in spl]
        return '\n'.join(spl)

    def __init__(self, mover=None, details=None):
        self._collapsed_samples = None
        self._samples = None
        self._len = None
        self._accepted = None
        self._accepted = None
        self.mover = mover
        self.subchanges = list()
        self.details = details

    def to_dict(self):
        return {
            'mover' : self.mover,
            'details' : self.details
        }

    @property
    def subchange(self):
        if len(self.subchanges) == 1:
            return self.subchanges[0]
        else:
            # TODO: might raise exception
            return None

    @property
    def opened(self):
        """
        Return the full PathMoveChange object

        Notes
        -----
        The full Movepath can only be returned if it has been generated and
        is still in memory. A collapsed Movepath that is loaded does NOT
        contain the full information anymore. This is the whole purpose of
        only storing essential information.
        """
        return self

    @property
    def closed(self):
        """
        Return a closed PathMoveChange object copy

        This will return a PathMoveChange object that will still know the used
        mover and all relevant samples. All underlying information about the
        move is hidden and will not be stored
        """
        obj = CollapsedPathMoveChange(trials=self.collapsed_samples, mover=self.mover)
        obj._subchange = self

        return obj

    def __iter__(self):
        yield self
        for subchange in self.subchanges:
            for change in subchange:
                yield change

    def __getitem__(self, item):
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
            return self._check_head_node(item)

            # Disable checking for submoves for now

            # the head node did not fit so continue trying subnodes
#            for sub in self.subchanges:
#                if item in sub:
#                    return True

            return False

        else:
            raise ValueError('Only PathMovers or PathMoveChanges can be tested.')

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

    def reduced(self, selected_samples = None, use_all_samples=False):
        """
        Reduce the underlying PathMoveChange to a subset of relevant samples

        This can be used to reduce a PathMoveChange (and everything below) to
        a subchange that only uses specific samples that can be picked
        from the list of all returned samples. Note that this can be
        problematic since a move might not return a different number of
        samples each time it is run.

        Parameters
        ----------
        selected_samples : list of int
            list of integer indices to be kept from the list of returned
            samples in this subchange. Slicing is not allowed, but negative
            indices are and conforms to the usual python convention (-1
            is the last samples, etc.)
        use_all_samples : bool
            if `True` the selected samples will be chosen from the list of
            all created samples (accepted and rejected ones), otherwise only
            the accepted ones will be chosen from. This is the default and
            corresponds to chose from `.samples`
        """

        # @DWHS do we want to allow also to select rejected samples? Meaning,
        # that we could allow the user to pick samples from .all_samples
        # or .samples

        # the check for collapsed_samples makes sure that at least the
        # relevant ones are present

        if selected_samples is None:
            # this case corresponds to .closed
            return self.closed

        if use_all_samples:
            # choose all generated samples
            sample_set = self.all_samples
        else:
            # chose only accepted ones!
            sample_set = self.samples

        # allow for negative indices to be picked, e.g. -1 is the last sample
        if len(selected_samples) > 0:
            selected_samples = [ idx % len(sample_set) if idx < 0 else idx for idx in selected_samples]

        samples = [
            samp for idx, samp in enumerate(sample_set)
            if idx in selected_samples or samp in self.collapsed_samples
        ]

        obj = CollapsedPathMoveChange(trials=samples, mover=self.mover)
        obj._subchange = self

        return obj

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
        if self._collapsed_samples is None:
            s = paths.SampleSet([]).apply_samples(self.samples)

            # keep order just for being thorough
            self._collapsed_samples = [
                samp for samp in self.samples
                if samp in s
            ]

        return self._collapsed_samples

    @property
    def trials(self):
        """
        A list of the samples that are needed to update the new sample_set

        Notes
        -----
        These samples are purely for this PathMoveChange node and not for any
        underlying nodes. It is effectively only used by SamplePathMoveChange
        """
        return self._trials

    @property
    def all_samples(self):
        """
        Returns a list of all samples generated during the PathMove.

        This includes all accepted and rejected samples (which does NOT
        include hidden samples yet)

        TODO: Decide, if we want hidden samples to be included or not
        TODO: might be obsolete
        """
        return self._trials

    @property
    def accepted(self):
        """
        Returns if this particular move was accepted.

        Mainly used for rejected samples.
        """
        if self._accepted is None:
            self._accepted = len(self.samples) > 0

        return self._accepted

    def __add__(self, other):
        return SequentialPathMoveChange([self, other])

    def apply_to(self, other):
        """
        Standard apply is to apply the list of samples contained
        """
        return other

    def __call__(self, other):
        return self.apply_to(other)

    def getkey(self, key):
        if self.mover is key[0]:
            subid = 0
            ret = self
            if len(key) > 1:
                for submove in zip(key[1][0::2], key[1][1::2]):
                    if submove[0] is not None:
                        subres = self.subchanges[subid]
                        if subres is not None:
                            ret = subres
                        else:
                            return None
            return ret
        else:
            return None

    @property
    def samples(self):
        """
        Returns a list of all samples that are accepted in this move

        This contains unnecessary, but accepted samples, too.
        """
        if self._samples is None:
            self._samples = self._get_samples()

        return self._samples

    def _get_samples(self):
        """
        Determines all relevant accepted samples for this move

        Includes all accepted samples also from submoves
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

    def apply_to(self, other):
        return other

    def __str__(self):
        return ''

    def to_dict(self):
        return {}

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
    def __init__(self, trials, mover=None, details=None):
        super(SamplePathMoveChange, self).__init__(mover=mover, details=details)
        self._trials = []

        if type(trials) is paths.Sample:
            trials = [trials]

        self._trials.extend(trials)

    def _get_samples(self):
        return [ sample for sample in self.trials if sample.accepted ]


    def to_dict(self):
        return {
            'mover' : self.mover,
            'trials' : self.trials,
            'details' : self.details
        }

    def apply_to(self, other):
        """
        Standard apply is to apply the list of samples contained
        """
        new_sample_set = paths.SampleSet(other).apply_samples(self.samples)
        #TODO: add pathmove change to new sampleset?
#        new_sample_set.pathmovechange = self
        return new_sample_set


@ops_object
class CollapsedPathMoveChange(SamplePathMoveChange):
    """
    Represent a collapsed PathMoveChange that has potential hidden sub moves
    """
    @property
    def opened(self):
        if hasattr(self, '_subchange') and self._subchange is not None:
            return self._subchange
        else:
            return self

    def closed(self):
        return self

    @property
    def collapsed_samples(self):
        if self._collapsed_samples is None:
            self._collapsed_samples = self.trials

        return self._collapsed_samples

    def __str__(self):
        if self.mover is not None:
            return '%s [collapsed] : %d samples' % (self.mover.cls, len(self.trials)) + ' ' + str(self.trials) + ''
        else:
            return '%s [collapsed] : %d samples' % ('CollapsedMove', len(self.trials)) + ' ' + str(self.trials) + ''


@ops_object
class RandomChoicePathMoveChange(PathMoveChange):
    """
    A PathMoveChange that represents the application of a mover chosen randomly
    """
    def __init__(self, subchange, mover=None, details=None):
        super(RandomChoicePathMoveChange, self).__init__(mover=mover, details=details)
        self.subchanges = [subchange]

    def to_dict(self):
        return {
            'mover' : self.mover,
            'subchange' : self.subchanges[0],
            'details' : self.details
        }

    def _get_samples(self):
        return self.subchange.samples

    def apply_to(self, other):
        return self.subchange.apply_to(other)

    def __str__(self):
        return 'RandomChoice :\n' + PathMoveChange._indent(str(self.subchange))

    @property
    def all_samples(self):
        return self.subchange.all_samples


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

    def _get_samples(self):
        samples = []
        for subchange in self.subchanges:
            samples = samples + subchange.samples
        return samples

    @property
    def all_samples(self):
        samples = []
        for subchange in self.subchanges:
            samples = samples + subchange.all_samples
        return samples

    def apply_to(self, other):
        sample_set = other

        for subchange in self.subchanges:
            sample_set = subchange.apply_to(sample_set)

        return sample_set

    def __str__(self):
        return 'SequentialMove : %s : %d samples\n' % \
               (self.accepted, len(self.samples)) + \
               PathMoveChange._indent('\n'.join(map(str, self.subchanges)))

@ops_object
class PartialAcceptanceSequentialPathMoveChange(SequentialPathMoveChange):
    """
    PartialAcceptanceSequentialMovePath has no own samples, only inferred
    Sampled from the underlying MovePaths
    """

    def _get_samples(self):
        changes = []
        for subchange in self.subchanges:
            if subchange.accepted:
                changes.extend(subchange.samples)
            else:
                break

        return changes

    #TODO: Could be removed
    def apply_to(self, other):
        sample_set = other

        for subchange in self.subchanges:
            if subchange.accepted:
                sample_set = subchange.apply_to(sample_set)
            else:
                break

        return sample_set

    def __str__(self):
        return 'PartialAcceptanceMove : %s : %d samples\n' % \
               (self.accepted, len(self.samples)) + \
               PathMoveChange._indent('\n'.join(map(str, self.subchanges)))

@ops_object
class ConditionalSequentialPathMoveChange(SequentialPathMoveChange):
    """
    ConditionalSequentialMovePath has no own samples, only inferred Samples
    from the underlying MovePaths
    """

    def apply_to(self, other):
        sample_set = other

        for subchange in self.subchanges:
            if subchange.accepted:
                sample_set = subchange.apply_to(sample_set)
            else:
                return other

        return sample_set

    def _get_samples(self):
        changes = []
        for subchange in self.subchanges:
            if subchange.accepted:
                changes.extend(subchange.samples)
            else:
                return []

        return changes

    def __str__(self):
        return 'ConditionalSequentialMove : %s : %d samples\n' % \
               (self.accepted, len(self.samples)) + \
               PathMoveChange._indent( '\n'.join(map(str, self.subchanges)))

@ops_object
class FilterSamplesPathMoveChange(PathMoveChange):
    """
    A PathMoveChange that keeps a selection of the underlying samples
    """

    def __init__(self, subchange, selected_samples, use_all_samples=False, mover=None, details=None):
        super(FilterSamplesPathMoveChange, self).__init__(mover=mover, details=details)
        self.subchanges = [subchange]
        self.selected_samples = selected_samples
        self.use_all_samples = use_all_samples


    def _get_samples(self):
        if self.use_all_samples:
            # choose all generated samples
            sample_set = self.all_samples
        else:
            # chose only accepted ones!
            sample_set = self.samples

        # allow for negative indices to be picked, e.g. -1 is the last sample
        samples = [ idx % len(sample_set) for idx in self.selected_samples]

        return samples

    def __str__(self):
        return 'FilterMove : pick samples [%s] from sub moves : %s : %d samples\n' % \
               (str(self.selected_samples), self.accepted, len(self.samples)) + \
               PathMoveChange._indent( str(self.subchange) )

    def to_dict(self):
        return {
            'subchange' : self.subchanges[0],
            'selected_samples' : self.selected_samples,
            'use_all_samples' : self.use_all_samples,
            'mover' : self.mover,
            'details' : self.details
        }

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

    def _get_samples(self):
        samples = self.subchange.samples
        if len(samples) > 1:
            samples = [samples[-1]]

        return samples

    def __str__(self):
        return 'Restrict to last sample : %s : %d samples\n' % \
               (self.accepted, len(self.samples)) + \
               PathMoveChange._indent( str(self.subchange) )

    def to_dict(self):
        return {
            'subchange' : self.subchanges[0],
            'mover' : self.mover,
            'details' : self.details
        }

@ops_object
class PathSimulatorPathMoveChange(PathMoveChange):
    """
    A PathMoveChange that just wraps a subchange and references a PathSimulator
    """

    def __init__(self, subchange, pathsimulator=None, step=-1, mover=None, details=None):
        super(PathSimulatorPathMoveChange, self).__init__(mover=mover, details=details)
        self.subchanges = [subchange]
        self.pathsimulator = pathsimulator
        self.step = step

    def _get_samples(self):
        return self.subchange.samples

    def __str__(self):
        return 'PathSimulatorStep : %s : Step # %d with %d samples\n' % \
               (str(self.pathsimulator.cls), self.step, len(self.samples)) + \
               PathMoveChange._indent( str(self.subchange) )

    def to_dict(self):
        return {
            'subchange' : self.subchanges[0],
            'pathsimulator' : self.pathsimulator,
            'step' : self.step,
            'mover' : self.mover,
            'details' : self.details
        }
