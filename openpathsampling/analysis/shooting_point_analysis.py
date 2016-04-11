import openpathsampling as paths

# based on http://stackoverflow.com/a/3387975
import collections
class TransformedDict(collections.MutableMapping):
    """A dictionary that applies an arbitrary key-altering
       function before accessing the keys"""

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
        return TransformedDict(new_hash, 
                               {self.hash_representatives[k]: self.store[k] 
                                for k in self.store})


class SnapshotByCoordinateDict(TransformedDict):
    def __init__(self, *args, **kwargs):
        hash_fcn = lambda x : x.coordinates.tostring()
        super(TransformedDict, self).__init__(hash_fcn, args, kwargs)


def shooting_point_analysis(steps, states):
    results = {}
    for step in steps:
        # TODO: this should in step.change.canonical.details
        details = step.change.canonical.trials[0].details
        try:
            shooting_snap = shooting_snapshot
        except AttributeError:
            # wrong kind of move
            pass
        except IndexError:
            # very wrong kind of move
            pass
        else:
            # easy to change how we define the key
            key = shooting_snap.coordinates.tostring()
            trial_traj = step.change.canonical.trials[0].trajectory
            init_traj = details.initial_trajectory
            shooting_traj = trial_traj.unique_subtrajectory(init_traj)
            endpoints = list(set([shooting_traj[0], shooting_snap[-1]]))
            # we use set in case there's only one frame (`first is last`)
            local = {}
            for state in states:
                winners = [snap for snap in endpoints if state(snap)]
                if len(winners) == 1:
                    # this is messy... there has to be a better way
                    try:
                        result_dict = results[key]
                    except KeyError:
                        results[key] = {state : 0 for state in states}
                        result_dict = results[key]
                    finally:
                        results[key][state] += 1
                elif len(winners) > 1:
                    print step.change
                    print winners
                    print trial_traj, init_traj, shooting_traj
                    raise RuntimeError
    
    return results

