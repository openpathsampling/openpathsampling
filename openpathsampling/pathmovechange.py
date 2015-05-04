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
        obj = CollapsedPathMoveChange(samples=self.collapsed_samples, mover=self.mover)
        obj._subchange = self

        return obj

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

    def _check_head_node(self, items):
        if isinstance(items[0], paths.PathMover):
            # a subtree of pathmovers
            if self.mover is items[0]:
                #print 'found head'
                # found current head node, check, if children match in order
                left = 0
                submovers = [ch.mover for ch in self.subchanges]
                subvalues = items[1]
                if type(subvalues) is not list:
                    subvalues = [subvalues]

                for sub in zip(subvalues[0::2], subvalues[1::2]):
                    if left >= len(self.subchanges):
                        # no more subchanges to match
                        return False
                    if sub is None:
                        # None is a placeholder so move token +1
                        left = left + 1
                    if type(sub) is dict:
                        if sub[0] is None:
                            while left < len(self.subchanges):
                                if not [self.subchanges[left].mover, [sub[1]]] in self.subchanges[left]:
                                    left = left + 1
                                else:
                                    left = left + 1
                                    break

                            if left == len(self.subchanges):
                                return False
                        elif sub[0] not in submovers[left:]:
                            #print 'missing sub', sub.keys()[0], 'in', submovers[left:]
                            return False
                        else:
                            idx = submovers.index(sub[0])
                            left = idx + 1
                            if not [sub[0], [sub[1]]] in self.subchanges[idx]:
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

        1. Subchanges are given using a dict { parent : child }
        2. Several subchanges are given in a list. [child1, child2]
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
#            for sub in self.subchanges:
#                if item in sub:
#                    return True

            return False

        else:
            raise ValueError('Only PathMovers or PathMoveChanges can be tested.')

    def tree(self):
        return {self : [ ch.tree() for ch in self.subchanges] }

    def movetree(self):
        return {self.mover : [ ch.movetree() for ch in self.subchanges] }

    def keytree(self, movepath=None):

        if movepath is None:
            movepath = [self.mover]

        result = list()
        result.append( ( movepath, self ) )
        mp = []
        for sub in self.subchanges:
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

        if len(self.subchanges) > 1:
            return { fnc(self, **kwargs) : [node.map_tree(fnc, **kwargs) for node in self.subchanges]}
        elif len(self.subchanges) == 1:
            return { fnc(self, **kwargs) : self.subchanges[0].map_tree(fnc, **kwargs)}
        else:
            return fnc(self, **kwargs)

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
            sample_set = self.results

        # allow for negative indices to be picked, e.g. -1 is the last sample
        if len(selected_samples) > 0:
            selected_samples = [ idx % len(sample_set) if idx < 0 else idx for idx in selected_samples]

        samples = [
            samp for idx, samp in enumerate(sample_set)
            if idx in selected_samples or samp in self.collapsed_samples
        ]

        obj = CollapsedPathMoveChange(samples=samples, mover=self.mover)
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
        if self._collapsed is None:
            self._collapsed = self.trials

        return self._collapsed

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