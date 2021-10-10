"""
This module includes tools to implement various common patterns used in
programming.
"""

import openpathsampling as paths

def listify(obj):
    """Force ``obj`` to be a list.
    """
    try:
        _ = iter(obj)
    except TypeError:
        obj = [obj]
        listified = True
    else:
        listified = False
    return obj, listified

def unlistify(iterable, listified):
    if listified:
        assert len(iterable) == 1, "Object to unlistify should be length 1"
        return iterable[0]
    else:
        return iterable


def trajectorify(thing):
    """
    Take a snapshot, trajectory, or iterable, and make a trajectory of it.

    This is intended to allow snapshots to be used in situations where we
    normally require trajectories. This does *not* provide information to
    "untrajectorify."
    """
    if isinstance(thing, paths.engines.Trajectory):
        iterable = thing
    elif hasattr(thing, '__iter__'):  # TODO: is this just isiterable?
        iterable = thing
    elif isinstance(thing, paths.engines.BaseSnapshot):
        iterable = [thing]
    else:
        raise RuntimeError("Can't make a trajectory out of " + str(thing))

    return paths.Trajectory(iterable)
