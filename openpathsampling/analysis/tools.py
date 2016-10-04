import logging

logger = logging.getLogger(__name__)

def pathlength(sample):
    return len(sample.trajectory)

def max_lambdas(sample, orderparameter):
    return max(orderparameter(sample.trajectory))

def sampleset_sample_generator(steps):
    for step in steps:
        sset = step.active # take the sampleset after the move
        for sample in sset:
            yield sample

def guess_interface_lambda(crossing_probability, direction=1):
    """
    Guesses the lambda for the interface based on the crossing probability.
    The assumption is that the interface lambda value is the last value
    where the (reverse) cumulative crossing probability is (nearly) 1.

    Parameters
    ----------
    crossing_probability : Histogram
        the max_lambda histogram
    direction : int
        if direction > 0, the order parameter is increasing, and the reverse
        cumulative histogram is used for the crossing probability. If
        direction < 0, the cumulative histogram is used for the crossing
        probability.

    Returns
    -------
    float
        the value of lambda for the interface
    """
    lambda_bin = -1
    if direction > 0:
        cp_vals = crossing_probability.reverse_cumulative().values()
        while abs(cp_vals[lambda_bin+1] - 1.0) < 1e-10:
            lambda_bin += 1
        outer_lambda = crossing_probability.bins[lambda_bin]
    elif direction < 0:
        cp_vals = crossing_probability.cumulative().values()
        while abs(cp_vals[lambda_bin+1] - 1.0) > 1e-10:
            lambda_bin += 1
        outer_lambda = crossing_probability.bins[lambda_bin-1]
    else:
        raise RuntimeError("Bad direction in guess_interface_lambda: " +
                           repr(direction))
    return outer_lambda


def minus_sides_summary(trajectory, minus_ensemble):
    # note: while this could be refactored so vol_dict is external, I don't
    # think this hurts speed very much, and it it really useful for testing
    minus_state = minus_ensemble.state_vol
    minus_innermost = minus_ensemble.innermost_vol
    minus_interstitial = minus_innermost & ~minus_state
    vol_dict = {
        "A" : minus_state,
        "X" : ~minus_innermost,
        "I" : minus_interstitial
    }
    summary = trajectory.summarize_by_volumes(vol_dict)
    # this is a per-trajectory loop
    count_sides = {"in" : [], "out" : []}
    side = None
    local_count = 0
    # strip off the beginning and ending in A
    for (label, count) in summary[1:-1]:
        if label == "X" and side != "out":
            if side == "in":
                count_sides["in"].append(local_count)
            side = "out"
            local_count = 0
        elif label == "A" and side != "in":
            if side == "out":
                count_sides["out"].append(local_count)
            side = "in"
            local_count = 0
        local_count += count
    return count_sides
