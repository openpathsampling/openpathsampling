import openpathsampling as paths
import collections
import pandas as pd

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

def shooting_point_analysis(steps, states):
    results = SnapshotByCoordinateDict()
    for step in steps:
        try:
            # TODO: this should in step.change.canonical.details
            details = step.change.canonical.trials[0].details
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
            trial_traj = step.change.canonical.trials[0].trajectory
            init_traj = details.initial_trajectory
            shooting_traj = trial_traj.unique_subtrajectory(init_traj)
            endpoints = list(set([shooting_traj[0], shooting_traj[-1]]))
            # we use set in case there's only one frame (`first is last`)
            initial = collections.Counter(
                {state: int(state(shooting_traj[0])) for state in states}
            )
            final = collections.Counter(
                {state: int(state(shooting_traj[-1])) for state in states}
            )
            total = initial + final if len(endpoints) == 2 else initial
            total_count = sum(total.values())
            assert total_count == 1 or total_count == 2
            try:
                results[key] += total
            except KeyError:
                results[key] = total

    return results

def shooting_point_analysis_to_pandas(results):
    """Each snapshot is a row, each state is a column"""
    transposed = pd.DataFrame(results.store).transpose().to_dict()
    df = pd.DataFrame(transposed)
    df.columns = [s.name for s in transposed.keys()]
    return df
