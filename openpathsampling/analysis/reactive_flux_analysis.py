import openpathsampling as paths
import collections
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from openpathsampling.analysis import ShootingPointAnalysis

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
        # ReactiveFluxAnalysis is inherited from ShootingPointAnalysis but
        # overrides the __init__ function. Thus, the call to the
        # SnapshotByCoordinateDict __init__ function has to be repeated
        # here (ugly)!
        paths.SnapshotByCoordinateDict.__init__(self)
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
            out_key = label_function(k)

            # look up the collection.Counter object for this snapshot
            counter_k = self[k]

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
        ndim = self._get_key_dim(list(r_store.keys())[0])

        # create histograms
        if ndim == 1:
            (shots_hist, b) = np.histogram(list(shots_dict.keys()),
                                           weights=list(shots_dict.values()),
                                           bins=bins)
            (sflux_hist, b) = np.histogram(list(sflux_dict.keys()),
                                           weights=list(sflux_dict.values()),
                                           bins=bins)
            b_list = [b]
        elif ndim == 2:
            (shots_hist, b_x, b_y) = np.histogram2d(
                x=[k[0] for k in shots_dict],
                y=[k[1] for k in shots_dict],
                weights=list(shots_dict.values()),
                bins=bins
            )
            (sflux_hist, b_x, b_y) = np.histogram2d(
                x=[k[0] for k in sflux_dict],
                y=[k[1] for k in sflux_dict],
                weights=list(sflux_dict.values()),
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
