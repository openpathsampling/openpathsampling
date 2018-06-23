import logging
import openpathsampling as paths
import collections
import numpy as np

from .shooting_point_analysis import SnapshotByCoordinateDict
from .shooting_point_analysis import ShootingPointAnalysis

logging.basicConfig()
logger = logging.getLogger(__name__)

class SShootingAnalysis(ShootingPointAnalysis):
    """
    Container and methods for S-shooting analysis.

    Parameters
    ----------
    steps : iterable of :class:`.MCStep` or None
        input MC steps to analyze; if None, no analysis performed
    states : list of :class:`.Volume`
        the volumes representing the stable states and the S region, in this
        order: [A, B, S].
    bias : :class:`.CollectiveVariable`
        bias weight function.
        
    """
    def __init__(self, steps, states, bias=None):
        # SShootingAnalysis is inherited from ShootinPointAnalysis but
        # overrides the __init__ function. Thus, the call to the
        # SnapshotByCoordinateDict __init__ function has to be repeated
        # here (ugly)!
        SnapshotByCoordinateDict.__init__(self)
        #super(ShootingPointAnalysis, self).__init__()

        # Definition of states A, B and S.
        self.state_A = states[0]
        self.state_B = states[1]
        self.state_S = states[2]

        # If called without bias create a dummy one.
        if bias is None:
            def dummy_bias(x):
                return 1.0
            cv_b = paths.CoordinateFunctionCV(name="cv_b", f=dummy_bias)
            self.bias = cv_b
        else:
            self.bias = bias

        # Initialize inconsistency counter.
        self.count_inconsistent_S = 0

        # Set one-direction shot length l.
        self.l = steps[0].simulation.trajectory_length
        
        # Prepare array for time correlation function.
        self._CABxhA_hS = None

        # Analyze each trajectory.
        if steps is not None:
            self.analyze(steps)

        # Issue a warning if some trajectories were harvested with a different
        # settings for the S region.
        if self.count_inconsistent_S > 0:
            logging.warning("Some trajectories ("
                            + str(self.count_inconsistent_S) +
                            " cases) were harvested with a different S region "
                            "definition than the one provided for the "
                            "analysis. Trajectories with starting snapshots "
                            "out of the given S region were ignored.")

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

        # Get key for dictionary entry for this starting snapshot.
        key = self.step_key(step)

        if key is not None:
            # Shortcut for current full trajectory (2*l+1 snapshots).
            trajectory = step.change.trials[-1]

            # Shortcut for forward/backward shot length.
            l = self.l

            # Check if S region from simulation matches with the one provided as
            # an argument.
            if self.state_S != step.simulation.state_S:
                self.count_inconsistent_S += 1
                if not self.state_S(trajectory[l]):
                    return;

            # If dictionary is empty for this snapshot create entry with
            # defaults:
            # M ......... Counter for trajectory segments.
            # I_Bt ...... Sum of 1/Bt[x(tau)] for all trajectory segments.
            # Ns_Bt ..... Sum of Ns[x(tau)]/Bt[x(tau)] for all trajectory
            #             segments.
            # hAhB_Bt ... Array with hA(0)*hB(t)/Bt[x(tau)] for all times t and
            #             all segments.
            if not (key in self):
                self[key] = {"M" : 0,
                             "Ns" : 0.0,
                             "I_Bt" : 0.0,
                             "Ns_Bt" : 0.0,
                             "hAhB_Bt" : np.zeros(l + 1)}

            # Calculate counter for number of states in S (Ns) and bias sum (Bt)
            # for first subtrajectory.
            Ns = 0
            Bt = 0.0
            for snapshot in trajectory[0:l+1]:
                if self.state_S(snapshot):
                    Ns += 1
                    Bt += self.bias(snapshot)

            # Loop over all starting points = subtrajectories.
            for istart in range(l + 1):

                # Update counter for number of states in S (Ns) and bias sum
                # (Bt)
                if istart > 0:
                    # Remove values for lost frame (most left one of the
                    # previous trajectory).
                    if self.state_S(trajectory[istart - 1]):
                        Ns -= 1
                        Bt -= self.bias(trajectory[istart - 1])
                    # Add values for new gained frame (frame right of the
                    # previous trajectory).
                    if self.state_S(trajectory[istart + l]):
                        Ns += 1
                        Bt += self.bias(trajectory[istart + l])

                # Increment trajectory counter and B-average of Ns / Bt.
                self[key]["M"] += 1
                self[key]["Ns"] += float(Ns)
                self[key]["I_Bt"] += 1.0 / Bt
                self[key]["Ns_Bt"] += Ns / Bt

                # Loop over states in subtrajectory.
                if self.state_A(trajectory[istart]):
                    for istep in range(l + 1):
                        if self.state_B(trajectory[istart + istep]):
                            self[key]["hAhB_Bt"][istep] += 1.0 / Bt
                #            self[key]["hAhB_Bt"][istep] += 0.5 / Bt
                #elif self.state_B(trajectory[istart]):
                #    for istep in range(l + 1):
                #        if self.state_A(trajectory[istart + istep]):
                #            self[key]["hAhB_Bt"][istep] += 0.5 / Bt
            return [k for k in self[key].keys()]
        else:
            return {}
                             
    def calculate_averages(self, label_function=None):
        """Calculate point-by-point and globally averaged results.

        This calculates the average time correlation function <hA(0)hB(t)>
        and other useful quantities globally and per initial snapshot.

        Parameters
        ----------
        label_function : callable
            the keys for the dictionary that is returned are
            `label_function(snapshot)`; default `None` gives the snapshot as
            key.

        Returns
        -------
        tuple
            format is (M, Ns, I_Bt, Ns_Bt, hAhB_Bt, dict)
            M : float
                Total number of harvested trajectory segments.
            Ns : float
                Average number of trajectory points in S.
            I_Bt : float
                Average of inverse bias sum over all trajectories.
            Ns_Bt : float
                Average of N_S divided by bias sum over all trajectories.
            hAhB_Bt : array
                Array with average of hA(0)hB(t) divided by bias sum
                over all trajectories.
            dict : dictionary
                all the above results for each snapshot separately.
        """
        # Set default label_function to identity.
        if label_function is None:
            label_function = lambda s : s

        # Initialize result variables.
        M = 0
        Ns = 0.0
        I_Bt = 0.0
        Ns_Bt = 0.0
        hAhB_Bt = np.zeros(self.l + 1)
        results = {}

        # Loop over all keys (this is a dictionary with snapshots
        # represented by coordinate hashes!).
        for k in self:
            # Set the output key (if label_function=None it is the snapshot).
            out_key = label_function(self.hash_representatives[k])

            # Look up the dictionary with data collected for this snapshot.
            data_k = self.store[k]

            # Increment global results.
            M += data_k["M"]
            Ns += data_k["Ns"]
            I_Bt += data_k["I_Bt"]
            Ns_Bt += data_k["Ns_Bt"]
            hAhB_Bt += data_k["hAhB_Bt"]

            # Save per-snapshot results.
            results[out_key] = {"M" : data_k["M"],
                                "Ns" : data_k["Ns"] / data_k["M"],
                                "I_Bt" : data_k["I_Bt"] / data_k["M"],
                                "Ns_Bt" : data_k["Ns_Bt"] / data_k["M"],
                                "hAhB_Bt" : data_k["hAhB_Bt"] / data_k["M"]}
        # Average global results.
        Ns /= M
        I_Bt /= M
        Ns_Bt /= M
        hAhB_Bt /= M
        self._CABxhA_hS = hAhB_Bt / Ns_Bt * (self.l + 1)

        # Return global quantities and dictionary with per-snapshot values.
        return (M, Ns, I_Bt, Ns_Bt, hAhB_Bt, results)

    def C_AB(self, hS=1.0, hA=1.0):
        """Calculate time correlation function C_AB(t).

        Automatically calls calculate_averages() if results are not yet stored
        internally.

        Parameters
        ----------
        hS : float
            Equilibrium propability for the system being in S, e.g. <hS> from
            umbrella sampling.
        hA : float
            Equilibrium propability for the system being in A, e.g. <hA> from
            umbrella sampling.

        Returns
        -------
        array
            Time correlation function C_AB(t).
        """
        return self.CABxhA_hS * hS / hA

    @property
    def CABxhA_hS(self):
        if self._CABxhA_hS is not None:
            return self._CABxhA_hS
        else:
            self.calculate_averages()
            return self._CABxhA_hS

    def committor(self, *args, **kwargs):
        raise NotImplementedError("The S-shooting analysis does not"
            + " calculate the committor.")

    def committor_histogram(self, *args, **kwargs):
        raise NotImplementedError("The S-shooting analysis does not"
            + "calculate the committor.")
