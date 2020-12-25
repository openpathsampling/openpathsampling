import collections
import pandas as pd
import numpy as np

from openpathsampling.progress import SimpleProgress

try:
    from collections import abc
except ImportError:
    import collections as abc


# based on http://stackoverflow.com/a/3387975
class TransformedDict(abc.MutableMapping):
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
        return iter(self.hash_representatives.values())

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
        hash_fcn = lambda x: x.coordinates.tobytes()
        super(SnapshotByCoordinateDict, self).__init__(hash_fcn,
                                                       *args, **kwargs)


class ShootingPointAnalysis(SimpleProgress, SnapshotByCoordinateDict):
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
        for step in self.progress(steps):
            self.analyze_single_step(step)

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
            total = collections.Counter({step[1]: 1})
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
            label_function = lambda s: s
        results = {}
        for k in self:
            out_key = label_function(k)
            counter_k = self[k]
            committor = float(counter_k[state]) / sum([counter_k[s]
                                                       for s in self.states])
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
        count_all = {k: sum(r_store[k].values()) for k in r_store}
        count_state = {k: r_store[k][state] for k in r_store}
        ndim = self._get_key_dim(list(r_store.keys())[0])
        if ndim == 1:
            (all_hist, b) = np.histogram(list(count_all.keys()),
                                         weights=list(count_all.values()),
                                         bins=bins)
            (state_hist, b) = np.histogram(list(count_state.keys()),
                                           weights=list(count_state.values()),
                                           bins=bins)
            b_list = [b]
        elif ndim == 2:
            (all_hist, b_x, b_y) = np.histogram2d(
                x=[k[0] for k in count_all],
                y=[k[1] for k in count_all],
                weights=list(count_all.values()),
                bins=bins
            )
            (state_hist, b_x, b_y) = np.histogram2d(
                x=[k[0] for k in count_state],
                y=[k[1] for k in count_state],
                weights=list(count_state.values()),
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
