import openpathsampling as paths
import collections
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# based on http://stackoverflow.com/a/3387975
class TransformedDict(collections.MutableMapping):
    """A dictionary that applies an arbitrary key-altering function before
    accessing the keys

    This implementation involves a particular hashing function. It is
    assumed that any two input objects which give the same hash are
    effectively identical, allowing later rehashing based on the same.
    """

    def __init__(self, hash_function, *args, **kwargs):
        self.store = dict()
        self.hash_representatives = dict()
        self.hash_function = hash_function
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        return self.store[self.hash_function(key)]

    def __setitem__(self, key, value):
        hashed = self.hash_function(key)
        if hashed not in self.hash_representatives:
            self.hash_representatives[hashed] = key
        self.store[hashed] = value

    def __delitem__(self, key):
        hashed = self.hash_function(key)
        del self.store[hashed]
        del self.hash_representatives[hashed]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def rehash(self, new_hash):
        """Create a new TransformedDict with this data and new hash.

        It is up to the user to ensure that the mapping from the old hash to
        the new is a function (i.e., each entry from the old hash can be
        mapped directly onto the new hash).

        For example, this is used to map from a snapshot's coordinates to
        a collective variable based on the coordinates. However, if the
        orignal hash was based on coordinates, but the new hash included
        velocities, the resulting mapping would be invalid. It is up to the
        user to avoid such invalid remappings.
        """
        return TransformedDict(new_hash, 
                               {self.hash_representatives[k]: self.store[k] 
                                for k in self.store})


class SnapshotByCoordinateDict(TransformedDict):
    """TransformedDict that uses snapshot coordinates as keys.

    This is primarily used to have a unique key for shooting point analysis
    (e.g., committor analysis).
    """
    def __init__(self, *args, **kwargs):
        hash_fcn = lambda x : x.coordinates.tostring()
        super(SnapshotByCoordinateDict, self).__init__(hash_fcn, 
                                                       *args, **kwargs)


class ShootingPointAnalysis(SnapshotByCoordinateDict):
    """
    Container and methods for shooting point analysis.

    This is especially useful for analyzing committors, which is
    automatically done on a per-configuration basis, and can also be done
    as a histogram.

    Parameters
    ----------
    steps : iterable of :class:`.MCStep` or None
        input MC steps to analyze; if None, no analysis performed
    states : list of :class:`.Volume`
        volumes to consider as states for the analysis. For pandas output,
        these volumes must be named.
    """
    def __init__(self, steps, states):
        super(ShootingPointAnalysis, self).__init__()
        self.states = states
        if steps is not None:
            self.analyze(steps)

    def analyze(self, steps):
        """Analyze a list of steps, adding to internal results.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep` or None
            MC steps to analyze
        """
        for step in steps:
            total = self.analyze_single_step(step)

    def analyze_single_step(self, step):
        """
        Analyzes final states from a path sampling step. Adds to internal
        results.

        Parameters
        ----------
        step : :class:`.MCStep`
            the step to analyze and add to this analysis

        Returns
        -------
        list of :class:`.Volume`
            the states which are identified as new final states from this
            move
        """
        key = self.step_key(step)
        if key is not None:
            details = step.change.canonical.details
            trial_traj = step.change.canonical.trials[0].trajectory
            init_traj = details.initial_trajectory
            test_points = [s for s in [trial_traj[0], trial_traj[-1]]
                           if s not in [init_traj[0], init_traj[-1]]]

            total = collections.Counter(
                {state: sum([int(state(pt)) for pt in test_points])
                            for state in self.states}
            )
            total_count = sum(total.values())
            # TODO: clarify assertion (at least one endpoint in state)
            assert total_count == 1 or total_count == 2
            try:
                self[key] += total
            except KeyError:
                self[key] = total
        else:
            total = {}

        return [s for s in total.keys() if total[s] > 0]

    @staticmethod
    def step_key(step):
        """
        Returns the key we use for hashing (the shooting snapshot).

        Parameters
        ----------
        step : :class:`.MCStep`
            the step to extract a shooting point from

        Returns
        -------
        :class:`.Snapshot` or None
            the shooting snapshot, or None if this step is not a shooting
            move.
        """
        key = None
        try:
            details = step.change.canonical.details
            shooting_snap = details.shooting_snapshot
        except AttributeError:
            # wrong kind of move (no shooting_snapshot)
            pass
        except IndexError:
            # very wrong kind of move (no trials!)
            pass
        else:
            # easy to change how we define the key
            key = shooting_snap
        return key

    @classmethod
    def from_individual_runs(cls, run_results, states=None):
        """Build shooting point analysis from pairs of shooting point to
        final state.

        Parameters
        ----------
        run_results : list of 2-tuples (:class:`.Snapshot`, :class:`.Volume`)
            the first element in each pair is the shooting point, the second
            is the final volume
        """
        if states is None:
            states = set(s[1] for s in run_results)
        analyzer = ShootingPointAnalysis(None, states)
        for step in run_results:
            key = step[0]
            total = collections.Counter({step[1] : 1})
            try:
                analyzer[key] += total
            except KeyError:
                analyzer[key] = total

        return analyzer

    def committor(self, state, label_function=None):
        """Calculate the (point-by-point) committor.

        This is for the point-by-point (per-configuration) committor, not
        for histograms. See `committor_histogram` for the histogram version.

        Parameters
        ----------
        state : :class:`.Volume`
            the committor is 1.0 if 100% of shots enter this state
        label_function : callable
            the keys for the dictionary that is returned are
            `label_function(snapshot)`; default `None` gives the snapshot as
            key.

        Returns
        -------
        dict :
            mapping labels given by label_function to the committor value
        """
        if label_function is None:
            label_function = lambda s : s
        results = {}
        for k in self:
            out_key = label_function(self.hash_representatives[k])
            counter_k = self.store[k]
            committor = float(counter_k[state]) / sum(counter_k.values())
            results[out_key] = committor
        return results

    @staticmethod
    def _get_key_dim(key):
        try:
            ndim = len(key)
        except TypeError:
            ndim = 1
        if ndim > 2 or ndim < 1:
            raise RuntimeError("Histogram key dimension {0} > 2 or {0} < 1 " 
                               + "(key: {1})".format(ndim, key))
        return ndim

    def committor_histogram(self, new_hash, state, bins=10):
        """Calculate the histogrammed version of the committor.

        Parameters
        ----------
        new_hash : callable
            values are histogrammed in bins based on new_hash(snapshot)
        state : :class:`.Volume`
            the committor is 1.0 if 100% of shots enter this state
        bins : see numpy.histogram
            bins input to numpy.histogram

        Returns
        -------
        tuple :
            hist, bins like numpy.histogram, where hist is the histogram
            count and bins is the bins output from numpy.histogram. 2-tuple
            in the case of 1D histogram, 3-tuple in the case of 2D histogram
        """
        rehashed = self.rehash(new_hash)
        r_store = rehashed.store
        count_all = {k : sum(r_store[k].values()) for k in r_store}
        count_state = {k : r_store[k][state] for k in r_store}
        ndim = self._get_key_dim(r_store.keys()[0])
        if ndim == 1:
            (all_hist, b) = np.histogram(count_all.keys(),
                                         weights=count_all.values(), 
                                         bins=bins)
            (state_hist, b) = np.histogram(count_state.keys(),
                                           weights=count_state.values(),
                                           bins=bins)
            b_list = [b]
        elif ndim == 2:
            (all_hist, b_x, b_y) = np.histogram2d(
                x=[k[0] for k in count_all],
                y=[k[1] for k in count_all],
                weights=count_all.values(),
                bins=bins
            )
            (state_hist, b_x, b_y) = np.histogram2d(
                x=[k[0] for k in count_state],
                y=[k[1] for k in count_state],
                weights=count_state.values(),
                bins=bins
            )
            b_list = [b_x, b_y]
        # if all_hist is 0, state_hist is NaN: ignore warning, return NaN
        with np.errstate(divide='ignore', invalid='ignore'):
            state_frac = np.true_divide(state_hist, all_hist)
        return tuple([state_frac] + b_list)

    def to_pandas(self, label_function=None):
        """
        Pandas dataframe. Row for each configuration, column for each state.

        Parameters
        ----------
        label_function : callable
            takes snapshot, returns index to use for pandas.DataFrame
        """
        transposed = pd.DataFrame(self.store).transpose().to_dict()
        df = pd.DataFrame(transposed)
        df.columns = [s.name for s in transposed.keys()]
        if label_function is None:
            df.index = range(len(df.index))
        else:
            # TODO: is ordering guaranteed here?
            df.index = [label_function(self.hash_representatives[k])
                        for k in self.store]
        return df


class ReactiveFluxAnalysis(ShootingPointAnalysis):
    """
    Container and methods for reactive flux analysis.

    Parameters
    ----------
    steps : iterable of :class:`.MCStep` or None
        input MC steps to analyze; if None, no analysis performed
    gradient : :class:`.CollectiveVariable`
        gradient of the reaction coordinate
    """

    def __init__(self, steps, gradient):
        # ReactiveFluxAnalysis is inherited from ShootinPointAnalysis but
        # overrides the __init__ function. Thus, the call to the
        # SnapshotByCoordinateDict __init__ function has to be repeated
        # here (ugly)!
        SnapshotByCoordinateDict.__init__(self)
        #super(ShootingPointAnalysis, self).__init__()
        self.gradient = gradient
        if steps is not None:
            self.analyze(steps)

    def analyze_single_step(self, step):
        """
        Adding a single step.

        Parameters
        ----------
        step : :class:`.MCStep`
            the step to analyze and add to this analysis

        Returns
        -------
        list of :class:`.Volume`
            the states which are identified as new final states from this
            move
        """
        key = self.step_key(step)
        if key is not None:
            # check if trajectory has been accepted
            accepted = step.change.canonical.accepted

            dldt = 0.0
            # if yes, calculate time derivative of reaction coordinate
            if accepted:
                # get shooting point and extract velocities
                snapshot = step.change.trials[0].trajectory[-1]
                v = snapshot.velocities

                # calculate gradient of reaction coordinate with respect
                # to shooting point coordinates
                dldx = self.gradient(snapshot)

                # calculate scalar product of gradient and velocities
                # = time derivative of reaction coordinate.
                # CAUTION: dldx and v are lists of vectors, for each particle
                # there is one pair of gradient and velocity vectors. Thus,
                # we have to sum over all particles in the end.
                dldt = sum([np.dot(a, b.T) for a, b in zip(dldx, v)])

            # set trial information in a collections.Counter dictionary
            total = collections.Counter(
                {"accepted" : int(accepted),
                 "rejected" : int(not accepted),
                 "sumflux"  : dldt}
            )
            # sum up Counter for this snapshot
            try:
                self[key].update(total)
            # if Counter does not exist create it
            except KeyError:
                self[key] = total
        else:
            total = {}

        return [s for s in total.keys()]

    def flux(self, label_function=None):
        """Calculate point-by-point and total flux.

        This calculates the average flux per initial snapshot and the total
        flux averaged over all of these configurations. Use `flux_histogram`
        for a histogram version.

        Parameters
        ----------
        label_function : callable
            the keys for the dictionary that is returned are
            `label_function(snapshot)`; default `None` gives the snapshot as
            key.

        Returns
        -------
        tuple
            format is (flux, dict)
            flux : float
                total flux averaged over all shooting points
            dict : dictionary
                flux per snapshot, mapping labels given by label_function
                to the flux value
        """
        # set default label_function to identity
        if label_function is None:
            label_function = lambda s : s

        # initialize result variables and counter
        results = {}
        flux = 0.0
        num_trials = 0

        # loop over all keys (this is a dictionary with snapshots
        # represented by coordinate hashes!)
        for k in self:
            # exclude the gradient function
            #if k is self.gradient:
            #    continue

            # set the output key (if label_function=None it is the snapshot)
            out_key = label_function(self.hash_representatives[k])

            # look up the collection.Counter object for this snapshot
            counter_k = self.store[k]

            # sum up the total flux and count trial moves
            flux = flux + counter_k["sumflux"]
            num_trials = num_trials + \
                         (counter_k["accepted"] + counter_k["rejected"])

            # calculate the average flux for this snapshot
            flux_per_snap = float(counter_k["sumflux"]) \
                            / (counter_k["accepted"] + counter_k["rejected"])

            # add the flux per snapshot to the output dictionary
            results[out_key] = flux_per_snap

        # return tuple with total flux and per-snapshot flux dicationary
        return (flux / float(num_trials), results)

    @staticmethod
    def _get_key_dim(key):
        try:
            ndim = len(key)
        except TypeError:
            ndim = 1
        if ndim > 2 or ndim < 1:
            raise RuntimeError("Histogram key dimension {0} > 2 or {0} < 1 "
                               + "(key: {1})".format(ndim, key))
        return ndim

    def flux_histogram(self, new_hash, bins=10):
        """Calculate the histogrammed version of the flux.

        Parameters
        ----------
        new_hash : callable
            values are histogrammed in bins based on new_hash(snapshot)
        bins : see numpy.histogram
            bins input to numpy.histogram

        Returns
        -------
        tuple :
            hist, bins like numpy.histogram, where hist is the average flux
            per histogram bin and bins is the bins output from numpy.histogram.
            2-tuple in the case of 1D histogram, 3-tuple in the case of 2D
            histogram
        """
        # rehash configurations according to provided hash function
        rehashed = self.rehash(new_hash)
        r_store = rehashed.store

        # collect total number of shots per configuration
        shots_dict = {k : r_store[k]["accepted"]
                        + r_store[k]["rejected"] for k in r_store}

        # collect total flux per configuration
        sflux_dict = {k : r_store[k]["sumflux"] for k in r_store}

        # determine dimension of histogram
        ndim = self._get_key_dim(r_store.keys()[0])

        # create histograms
        if ndim == 1:
            (shots_hist, b) = np.histogram(shots_dict.keys(),
                                           weights=shots_dict.values(),
                                           bins=bins)
            (sflux_hist, b) = np.histogram(sflux_dict.keys(),
                                           weights=sflux_dict.values(),
                                           bins=bins)
            b_list = [b]
        elif ndim == 2:
            (shots_hist, b_x, b_y) = np.histogram2d(
                x=[k[0] for k in shots_dict],
                y=[k[1] for k in shots_dict],
                weights=shots_dict.values(),
                bins=bins
            )
            (sflux_hist, b_x, b_y) = np.histogram2d(
                x=[k[0] for k in sflux_dict],
                y=[k[1] for k in sflux_dict],
                weights=sflux_dict.values(),
                bins=bins
            )
            b_list = [b_x, b_y]

        # divide histograms to determine final result
        flux_hist = np.ma.true_divide(
            sflux_hist,
            np.ma.array(shots_hist,mask=shots_hist==0)
        )

        return tuple([flux_hist] + b_list)

    def committor(self, *args, **kwargs):
        raise NotImplementedError("The reactive flux analysis does not"
            + " calculate the committor.")

    def committor_histogram(self, *args, **kwargs):
        raise NotImplementedError("The reactive flux analysis does not"
            + "calculate the committor.")
