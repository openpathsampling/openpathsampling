import random
import logging

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableObject, lazy_loading_attributes
from openpathsampling.netcdfplus import DelayedLoader

from openpathsampling.tools import refresh_output

from collections import Counter

import sys
if sys.version_info > (3, ):
    from collections.abc import Mapping
else:
    from collections import Mapping


logger = logging.getLogger(__name__)


class SampleKeyError(Exception):
    def __init__(self, key, sample, sample_key):
        self.key = key
        self.sample = sample
        self.sample_key = sample_key
        self.msg = (str(self.key) + " does not match " + str(self.sample_key)
                    + " from " + str(self.sample))


# @lazy_loading_attributes('movepath')
class SampleSet(StorableObject, Mapping):
    """
    SampleSet is essentially a list of samples, with a few conveniences.  It
    can be treated as a list of samples (using, e.g., .append), or as a
    dictionary of ensembles mapping to a list of samples, or as a dictionary
    of replica IDs to samples. Replica ID has to an integer but it can be
    negative or zero.

    The dictionaries ensemble_dict and replica_dict are conveniences which
    should be kept consistent by any method which modifies the container.
    They do not need to be stored.

    Note
    ----
        Current implementation is as an unordered set. Therefore we don't
        have some of the convenient tools in Python sequences (e.g.,
        slices). On the other hand, I'm not sure whether that is meaningful
        here.
        Since replicas are integers we add slicing/ranges for replicas. In
        addition we support any iterable as input in __getitem__ an it will
        return an iterable over the results. This makes it possible to write
        `sset[0:5]` to get a list of of ordered samples by replica_id, or
        sset[list_of_ensembles].  replica_ids can be any number do not have
        to be subsequent to slicing does not make sense and we ignore it. We
        will also ignore missing replica_ids. A slice `1:5` will return all
        existing replica ids >=1 and <5. If you want exactly all replicas
        from 1 to 4 use `sset[xrange(1,5)]`


    Attributes
    ----------
    samples : list of Sample
        The samples included in this set.
    ensemble_dict : dict
        A dictionary with Ensemble objects as keys and lists of Samples as
        values.
    replica_dict : dict
        A dictionary with replica IDs as keys and lists of Samples as values
    """

    movepath = DelayedLoader()

    def __init__(self, samples, movepath=None):
        super(SampleSet, self).__init__()

        self._lazy = {}

        self.samples = []
        self.ensemble_dict = {}
        self.replica_dict = {}
        self.extend(samples)
        self.movepath = movepath

    @property
    def ensembles(self):
        return self.ensemble_dict.keys()

    def values(self):
        return self.samples

    @property
    def replicas(self):
        return self.replica_dict.keys()

    def __getitem__(self, key):
        if isinstance(key, paths.Ensemble):
            return random.choice(self.ensemble_dict[key])
        elif type(key) is int:
            return random.choice(self.replica_dict[key])
        elif hasattr(key, '__iter__'):
            return (self[element] for element in key)
        elif type(key) is slice:
            rep_idxs = filter(
                lambda x:
                    (key.start is None or x >= key.start) and
                    (key.stop is None or x < key.stop),
                sorted(self.replica_dict.keys())

            )

            return (self[element] for element in rep_idxs)

    def __setitem__(self, key, value):
        # first, we check whether the key matches the sample: if no, KeyError
        if isinstance(key, paths.Ensemble):
            if key != value.ensemble:
                raise SampleKeyError(key, value, value.ensemble)
        else:
            if key != value.replica:
                raise SampleKeyError(key, value, value.replica)

        if value in self.samples:
            # if value is already in this, we don't need to do anything
            return
        # Setting works by replacing one with the same key. We pick one with
        # this key at random (using __getitem__), delete it, and then append
        # the new guy. If nothing exists with the desired key, this is the
        # same as append.
        try:
            dead_to_me = self[key]
        except KeyError:
            dead_to_me = None
        if dead_to_me is not None:
            del self[dead_to_me]
        self.append(value)

    def __eq__(self, other):
        return Counter(self.samples) == Counter(other.samples)

    def __ne__(self, other):
        return not self == other

    def __delitem__(self, sample):
        self.ensemble_dict[sample.ensemble].remove(sample)
        self.replica_dict[sample.replica].remove(sample)
        if len(self.ensemble_dict[sample.ensemble]) == 0:
            del self.ensemble_dict[sample.ensemble]
        if len(self.replica_dict[sample.replica]) == 0:
            del self.replica_dict[sample.replica]
        self.samples.remove(sample)

    # TODO: add support for remove and pop

    def __iter__(self):
        for sample in self.samples:
            yield sample

    def __len__(self):
        return len(self.samples)

    def __contains__(self, item):
        # check for Sample, replica (int) and Ensemble, too
        if item in self.samples:
            return True
        elif item in self.ensemble_dict:
            return True
        elif item in self.replica_dict:
            return True
        else:
            return False

    def all_from_ensemble(self, ensemble):
        try:
            return self.ensemble_dict[ensemble]
        except KeyError:
            return []

    def all_from_replica(self, replica):
        try:
            return self.replica_dict[replica]
        except KeyError:
            return []

    def append(self, sample):
        if sample in self.samples:
            # question: would it make sense to raise an error here? can't
            # have more than one copy of the same sample, but should we
            # ignore it silently or complain?
            return

        self.samples.append(sample)
        try:
            self.ensemble_dict[sample.ensemble].append(sample)
        except KeyError:
            self.ensemble_dict[sample.ensemble] = [sample]
        try:
            self.replica_dict[sample.replica].append(sample)
        except KeyError:
            self.replica_dict[sample.replica] = [sample]

    def extend(self, samples):
        # note that this works whether the parameter samples is a list of
        # samples or a SampleSet!
        if type(samples) is not paths.Sample and hasattr(samples, '__iter__'):
            for sample in samples:
                self.append(sample)
        else:
            # also acts as .append() if given a single sample
            self.append(samples)

    def apply_samples(self, samples, copy=True):
        """Update by setting samples by replica in the order given

        """
        if isinstance(samples, Sample):
            samples = [samples]
        elif isinstance(samples, paths.MoveChange):
            samples = samples.results
        if copy:
            newset = SampleSet(self)
        else:
            newset = self
        for sample in samples:
            if type(sample) is not paths.Sample:
                raise ValueError(
                    'No SAMPLE! Type `%s` found.' % sample.__class__.__name__)
            # TODO: should time be a property of Sample or SampleSet?
            newset[sample.replica] = sample
        return newset

    def replica_list(self):
        """Returns the list of replicas IDs in this SampleSet

        """
        return self.replica_dict.keys()

    def ensemble_list(self):
        """Returns the list of ensembles in this SampleSet

        """
        return self.ensemble_dict.keys()

    def sanity_check(self):
        """Checks that the sample trajectories satisfy their ensembles

        """
        logger.info("Starting sanity check")
        for sample in self:
            logger.info("Checking sanity of " + repr(sample.ensemble) +
                        " with " + str(sample.trajectory))
            try:
                assert(sample.ensemble(sample.trajectory))
            except AssertionError as e:
                failmsg = ("Trajectory does not match ensemble for replica "
                           + str(sample.replica))
                if not e.args:
                    e.args = [failmsg]
                else:
                    arg0 = failmsg + e.args[0]
                    e.args = tuple([arg0] + list(e.args[1:]))
                raise  # re-raise last exception

    def consistency_check(self):
        """Check that all internal dictionaries are consistent

        This is mainly a sanity check for use in testing, but might be
        good to run (rarely) in the code until we're sure the tests cover
        all use cases.
        """

        # check that we have the same number of samples in everything
        nsamps_ens = 0
        for ens in self.ensemble_dict.keys():
            nsamps_ens += len(self.ensemble_dict[ens])
        nsamps_rep = 0
        for rep in self.replica_dict.keys():
            nsamps_rep += len(self.replica_dict[rep])
        nsamps = len(self.samples)
        assert nsamps == nsamps_ens, \
            "nsamps != nsamps_ens : %d != %d" % (nsamps, nsamps_ens)
        assert nsamps == nsamps_rep, \
            "nsamps != nsamps_rep : %d != %d" % (nsamps, nsamps_rep)

        # if we have the same number of samples, then we check that each
        # sample in samples is in each of the dictionaries
        for samp in self.samples:
            assert samp in self.ensemble_dict[samp.ensemble], \
                "Sample not in ensemble_dict! %r %r" % (
                    samp, self.ensemble_dict)
            assert samp in self.replica_dict[samp.replica], \
                "Sample not in replica_dict! %r %r" % (samp, self.replica_dict)

        # finally, check to be sure that there are no duplicates in
        # self.samples; this completes the consistency check
        for samp in self.samples:
            assert self.samples.count(samp) == 1, \
                    "More than one instance of %r!" % samp

    def append_as_new_replica(self, sample):
        """
        Adds the given sample to this SampleSet, with a new replica ID.

        The new replica ID is taken to be one greater than the highest
        previous replica ID.
        """
        if len(self) == 0:
            max_replica = -1
        else:
            max_replica = max([s.replica for s in self.samples])
        self.append(Sample(
            replica=max_replica + 1,
            trajectory=sample.trajectory,
            ensemble=sample.ensemble,
            bias=sample.bias,
            # details=sample.details,
            parent=sample.parent,
            mover=sample.mover
        ))

    @staticmethod
    def map_trajectory_to_ensembles(trajectory, ensembles):
        """Return SampleSet mapping one trajectory to all ensembles.

        One common approach to starting a simulation is to take a single
        transition trajectory (which satisfies all ensembles) and use it as
        the starting point for all ensembles.
        """
        return SampleSet([
            Sample.initial_sample(
                replica=ensembles.index(e),
                trajectory=paths.Trajectory(trajectory.as_proxies()),  # copy
                ensemble=e)
            for e in ensembles
        ])

    @staticmethod
    def translate_ensembles(sset, new_ensembles):
        """Return SampleSet using `new_ensembles` as ensembles.

        This creates a SampleSet which replaces the ensembles in the old
        sampleset with equivalent ensembles from a given list. The string
        description of the ensemble is used as a test.

        Note that this assumes that there are no one-to-many or many-to-one
        relations in the ensembles. If there are, then there is no unique
        way to translate.

        The approach used here will return the SampleSet with the maximum
        number of ensembles that overlap between the two groups.
        """
        translation = {}
        for ens1 in sset.ensemble_list():
            for ens2 in new_ensembles:
                if ens1.__str__() == ens2.__str__():
                    translation[ens1] = ens2

        new_samples = []
        for ens in translation:
            old_samples = sset.all_from_ensemble(ens)
            for s in old_samples:
                new_samples.append(Sample(
                    replica=s.replica,
                    ensemble=translation[s.ensemble],
                    trajectory=s.trajectory
                ))
        res = SampleSet.relabel_replicas_per_ensemble(SampleSet(new_samples))
        return res

    @staticmethod
    def relabel_replicas_per_ensemble(ssets):
        """
        Return a SampleSet with one trajectory ID per ensemble in `ssets`

        This is used if you create several sample sets (e.g., from
        bootstrapping different transitions) which have the same trajectory
        ID associated with different ensembles.
        """
        if type(ssets) is SampleSet:
            ssets = [ssets]
        samples = []
        repid = 0
        for sset in ssets:
            for s in sset:
                samples.append(Sample(
                    replica=repid,
                    trajectory=s.trajectory,
                    ensemble=s.ensemble
                ))
                repid += 1
        return SampleSet(samples)

    def generate_from_trajectories(
            self,
            ensembles,
            trajectories,
            preconditions=None,
            strategies=None,
            reuse_strategy='avoid',
            engine=None):
        """
        Create a SampleSet with as many initial samples as possible.

        The goal of this is to give the initial SampleSet that would be
        desired.

        Parameters
        ----------
        trajectories : list of :class:`.Trajectory` or :class:`.Trajectory`
            the input trajectories to use
        ensembles : list of :class:`Ensemble` or list of :class:`Ensemble`
            the list of ensembles to be generated. If an element is itself a
            list then one sample for one of the ensembles in that list if
            generated
        preconditions : list of str
            a list of possible steps to modify the initial list of trajectories.
            possible choices are

                1.  `sort-shortest` - sorting by shortest first,
                2.  `sort_median` - sorting by the middle one first and then in
                    move away from the median length
                3.  `sort-longest` - sorting by the longest first
                4.  `reverse` - reverse the order and
                5.  `mirror` which will add the reversed trajectories to the
                    list in the same order

            Default is `None` which means to do nothing.
        strategies : dict
            a dict that specifies the options used when ensemble functions
            are used to create a new sample.
        reuse_strategy : str
            if `avoid` then in a second attempt the used trajectories are
            tried
        engine : :class:`openpathsampling.engines.DyanmicsEngine`
            the engine used for extending moves

        Returns
        -------
        :class:`.SampleSet`
            sampleset with samples for every initial ensemble for this
            scheme that could be satisfied by the given trajectories

        See Also
        --------
        list_initial_ensembles
        """

        implemented_strategies = [
            'get',              # look for existing trajectories
            'split',            # look for existing sub-trajectories
            'extend-complex',   # try to extend long sub-trajectories
            'extend-minimal'    # try to extend short sub-trajectories
        ]

        implemented_preconditions = [
            'sort-shortest',
            'sort-median',
            'sort-longest',
            'reverse',
            'mirror'
        ]

        if preconditions is None:
            preconditions = ['mirror']

        # create a list of trajectories
        trajectories = paths.Trajectory._to_list_of_trajectories(trajectories)

        for pre in preconditions:
            if pre not in implemented_preconditions:
                raise RuntimeError(
                    '%s is not a valid precondition strategy. Choose from %s.' %
                    (pre, implemented_preconditions)
                )

            if pre == 'sort-shortest':
                trajectories = sorted(trajectories, key=len)
            elif pre == 'sort-longest':
                trajectories = sorted(trajectories, key=len)
            elif pre == 'sort-median':
                sorted_trajectories = sorted(trajectories, key=len)
                trajectories = list([p for p2 in zip(
                    sorted_trajectories[len(sorted_trajectories) / 2:],
                    reversed(sorted_trajectories[:len(sorted_trajectories) / 2])
                ) for p in p2])

                if len(sorted_trajectories) & 1:
                    trajectories.append(sorted_trajectories[-1])
            elif pre == 'reverse':
                trajectories = list(reversed(trajectories))
            elif pre == 'mirror':
                trajectories = trajectories + \
                    [traj.reversed for traj in trajectories]

        # let's always try the short trajectories first
        # print map(lambda x: hex(id(x)), trajectories)

        # we will try forward/backward interleaved
        used_trajectories = []

        # if we start with an existing sample set look at what we got
        # if we avoid we move the used ones to the back of the list
        # if we remove we remove the used ones
        for s in self:
            traj = s.trajectory
            if traj in trajectories:
                used_trajectories.append(traj)

        used_trajectories = sorted(used_trajectories, key=len)

        # print map(lambda x: hex(id(x)), used_trajectories)

        ensembles = [[x] if type(x) is not list else x for x in ensembles]

        # 1. look in the existing sample_set
        ensembles_to_fill, extra_ensembles = self.check_ensembles(ensembles)

        # we reverse because we want to be able to remove elements
        # from the list as we discover samples. This is easier to do
        # when the list is traversed backwards since indices to not
        # change, hence we reverse the list of ensemble and then traverse it
        # in reversed order
        ensembles_to_fill = \
            list(reversed(ensembles_to_fill))

        # 2. try strategies
        if strategies is None:
            # this is the default
            strategies = [
                'get',
                'split'
            ]

        for idx, strategy in enumerate(strategies):
            if type(strategy) is str:
                strategies[idx] = (strategy, dict())

            if strategies[idx][0] not in implemented_strategies:
                raise RuntimeError(
                    'Strategy `%s` is not known. Chose from %s.' % (
                        strategies[idx][0],
                        implemented_strategies
                    )
                )

        found_samples_str = ''
        for pos, ens_list in enumerate(ensembles):
            found_samples_str += '.' if ens_list in ensembles_to_fill else '+'

        for str_idx, (strategy, options) in enumerate(strategies):
            for idx, ens_list in reversed(list(enumerate(ensembles_to_fill))):
                pos = ensembles.index(ens_list)

                found_samples_str = \
                    found_samples_str[:pos] + \
                    '?' + found_samples_str[pos + 1:]

                refresh_output((
                    '# trying strategy #%d `%s`: still missing %d samples\n'
                    '%s\n'
                ) % (
                        str_idx + 1,
                        strategy,
                        len(ensembles_to_fill),
                        found_samples_str
                    ), ipynb_display_only=True, print_anyway=False)
                if type(ens_list) is not list:
                    ens_list = [ens_list]

                found = False

                for ens in ens_list:
                    # create the list of options to be passed on
                    opts = {key: value for key, value in options.items()
                            if key not in ['exclude']}

                    # exclude contains the Ensemble classes to be ignored
                    if 'exclude' in options:
                        if isinstance(ens, options['exclude']):
                            continue

                    # fill only the first in ens_list that can be filled

                    sample = None

                    if strategy == 'get':
                        sample = ens.get_sample_from_trajectories(
                            trajectories=trajectories,
                            used_trajectories=used_trajectories,
                            reuse_strategy=reuse_strategy,
                            **opts
                        )
                    elif strategy == 'split':
                        sample = ens.split_sample_from_trajectories(
                            trajectories=trajectories,
                            used_trajectories=used_trajectories,
                            reuse_strategy=reuse_strategy,
                            **opts
                        )
                    elif strategy == 'extend-complex' and engine:
                        if hasattr(ens, 'extend_sample_from_trajectories'):
                            sample = ens.extend_sample_from_trajectories(
                                trajectories=trajectories,
                                engine=engine,
                                level='complex',
                                **opts
                            )
                    elif strategy == 'extend-minimal' and engine:
                        if hasattr(ens, 'extend_sample_from_trajectories'):
                            sample = ens.extend_sample_from_trajectories(
                                trajectories=trajectories,
                                engine=engine,
                                level='minimal',
                                **opts
                            )
                    elif strategy == 'extend-native' and engine:
                        if hasattr(ens, 'extend_sample_from_trajectories'):
                            sample = ens.extend_sample_from_trajectories(
                                trajectories=trajectories,
                                engine=engine,
                                level='native',
                                **opts
                            )

                    # now, if we've found a sample, add it and
                    # make sure we chose a proper replica ID
                    if sample is not None:
                        found = True

                        # another way would be to look for the smallest not
                        # taken id. This one is simpler
                        if len(self.replicas) > 0:
                            replica_idx = max(0, max(self.replicas) + 1)
                        else:
                            replica_idx = 0

                        sample.replica = replica_idx

                        logger.info((
                            'generating - ensemble `%s` found sample '
                            'replica %d, length %d\n')
                            % (
                                ens.name, sample.replica, len(sample)
                            ))

                        self.append(sample)
                        if reuse_strategy != 'all':
                            # we mark the trajectory and its reversed as used
                            if sample.trajectory not in used_trajectories and (
                                not reuse_strategy.endswith('symmetric') or
                                sample.trajectory.reversed in used_trajectories
                            ):
                                used_trajectories.append(sample.trajectory)
                            # if reuse_strategy.endswith('symmetric'):
                            #     used_trajectories.append(
                            #         sample.trajectory.reversed)

                            # we want the list of used_trajectories to be
                            # sorted. Short ones first. So if we have to chose
                            # from the used_ones, use the shortest one
                            # used_trajectories = sorted(
                            #     used_trajectories, key=len)

                        # found a sample in this category so remove it for
                        # other tries
                        del ensembles_to_fill[idx]

                        # do not try other ensembles in this category
                        break

                found_samples_str = \
                    found_samples_str[:pos] + \
                    (str(str_idx + 1)[0] if found else '.') + \
                    found_samples_str[pos + 1:]

        refresh_output((
            '# finished generating: still missing %d samples\n'
            '%s\n'
        ) % (
            len(ensembles_to_fill),
            found_samples_str
        ), ipynb_display_only=True, print_anyway=False)

        return self

    def check_ensembles(self, ensembles):
        """
        Check for missing or extra ensembles in the sampleset

        This is primary used programmatically as a reusable function for
        several use cases where we need this information. See functions
        under "see also" for examples of such cases.

        Parameters
        ----------
        ensembles : list of :class:`Ensemble` or list of :class:`Ensemble`
            the list of ensembles to be generated. If an element is itself a
            list then one sample for one of the ensembles in that list if
            generated

        Returns
        -------
        missing : list of list of :class:`.Ensemble`
            ensembles needed by the move scheme and missing in the sample
            set, in the format used by `list_initial_ensembles`
        extra : list of :class:`.Ensemble`
            ensembles in the sampleset that are not used by the

        See Also
        --------
        MoveScheme.list_initial_ensembles
        MoveScheme.assert_initial_conditions
        MoveScheme.initial_conditions_report
        """
        samples = paths.SampleSet(self)  # to make a copy
        missing = []
        for ens_list in ensembles:
            if type(ens_list) is not list:
                ens_list = [ens_list]
            sample = None
            for ens in ens_list:
                if ens in samples.ensemble_list():
                    sample = samples[ens]
                    break
            if sample is not None:
                del samples[sample]
            else:
                missing.append(ens_list)

        # missing, extra
        return missing, samples.ensemble_list()

    def copy_without_parents(self):
        """
        Return a copy of the sample set where all samples.parents are removed

        Useful, if you are not interested in the heritage of the sample and
        store a clean set of samples
        """
        return SampleSet(
            [Sample(
                replica=s.replica,
                trajectory=s.trajectory,
                ensemble=s.ensemble,
                bias=s.bias,
                # details=s.details,
                mover=s.mover
            ) for s in self]
        )


# @lazy_loading_attributes('parent', 'mover')
class Sample(StorableObject):
    """
    A Sample represents a given "draw" from its ensemble, and is the return
    object from a PathMover. It and contains all information about the move,
    initial trajectories, new trajectories (both as references).

    Since each Sample is a single representative of a single ensemble, each
    Sample consists of one replica ID, one trajectory, and one ensemble.
    This means that movers which generate more than one "draw" (often from
    different ensembles, e.g. replica exchange) will generate more than one
    Sample object.

    Attributes
    ----------
    replica : int
        The replica ID to which this Sample applies. The replica ID can also
        be negative.
    trajectory : openpathsampling.Trajectory
        The trajectory (path) for this sample
    ensemble : openpathsampling.Ensemble
        The Ensemble this sample is drawn from
    """

    parent = DelayedLoader()
    mover = DelayedLoader()

    def __init__(self,
                 replica=None,
                 trajectory=None,
                 ensemble=None,
                 bias=1.0,
                 parent=None,
                 mover=None
                 ):

        super(Sample, self).__init__()

        self._lazy = {}

        self.bias = bias
        self.replica = replica
        self.ensemble = ensemble
        self.trajectory = trajectory
        self.parent = parent
        # self.details = details
        self.mover = mover

        self._valid = None

    # ==========================================================================
    # LIST INHERITANCE FUNCTIONS
    # ==========================================================================

    def __len__(self):
        return len(self.trajectory)

    def __getslice__(self, *args, **kwargs):
        return self.trajectory.__getslice__(*args, **kwargs)

    def __getitem__(self, *args, **kwargs):
        return self.trajectory.__getitem__(*args, **kwargs)

    def __reversed__(self):
        """
        Return a reversed iterator over all snapshots in the samples trajectory

        Returns
        -------
        Iterator()
            The iterator that iterates the snapshots in reversed order

        Notes
        -----
        A reversed trajectory also has reversed snapshots! This means
        that Trajectory(list(reversed(traj))) will lead to a time-reversed
        trajectory not just frames in reversed order but also reversed momenta.

        """
        if self.trajectory is not None:
            return reversed(self.trajectory)
        else:
            return []  # empty iterator

    def as_proxies(self):
        if self.trajectory is not None:
            return self.trajectory.as_proxies()
        else:
            return []  # empty iterator

    def __iter__(self):
        """
        Return an iterator over all snapshots in the samples trajectory

        Returns
        -------
        Iterator()
            The iterator that iterates the snapshots

        """
        if self.trajectory is not None:
            return iter(self.trajectory)
        else:
            return []  # empty iterator

    def __str__(self):
        # mystr  = "Replica: "+str(self.replica)+"\n"
        # mystr += "Trajectory: "+str(self.trajectory)+"\n"
        # mystr += "Ensemble: "+repr(self.ensemble)+"\n"
        mystr = 'Sample(RepID: %d, Ens: %s, %s)' % (
            self.replica, repr(self.ensemble), repr(self.trajectory))
        return mystr

    @property
    def valid(self):
        """Returns true if a sample is in its ensemble

        Returns
        -------
        bool
            `True` if the trajectory is in the ensemble `False` otherwise
        """
        if self._valid is None:
            if self.trajectory is None:
                self._valid = True
            else:
                if self.ensemble is not None:
                    self._valid = self.ensemble(self.trajectory)
                else:
                    # no ensemble means ALL ???
                    self._valid = True

        return self._valid

    def __repr__(self):
        return '<Sample @ ' + str(hex(id(self))) + '>'

    def copy_reset(self):
        """
        Copy of Sample with initialization move details.
        """
        result = Sample(
            replica=self.replica,
            trajectory=self.trajectory,
            ensemble=self.ensemble
        )
        return result

    @staticmethod
    def initial_sample(replica, trajectory, ensemble):
        """
        Initial sample from scratch.

        Used to create sample in a given ensemble when generating initial
        conditions from trajectories.
        """
        result = Sample(
            replica=replica,
            trajectory=trajectory,
            ensemble=ensemble
        )
        return result

    @property
    def acceptance(self):
        if not self.valid:
            return 0.0

        return self.bias

    @property
    def heritage(self):
        samp = self
        while samp.parent is not None:
            # just one sample so use this
            yield samp
            samp = samp.parent
