#!/usr/bin/env python2

from openpathsampling.pathsimulator import PathSimulator
import openpathsampling as paths
import os


class TransitionStateEnsemble(PathSimulator):
    """
    Transition state ensemble simulation. Which snapshot from a trajectory
    results in a certain committor?

    Parameters
    ----------
    trajectories : list of :class:`.Trajectory`
        list of trajectories to get the transition state ensemble from
    stateA : :class:`.Volume`
        the volume representing the initial state of the transition
    stateB : :class:`.Volume`
        the volume representing the final state of the transition
    engine : :class:`.DynamicsEngine`
        the dynamics engine to use in the committor simulaion
    storage : :class:`.Storage`
        the file to store the simulations in
    randomizer : :class:`.SnapshotModifier`
        the method used to modify the input snapshot before each committor shot
    n_per_snapshot : int
        the amount of committor shots to run from each snapshot
    pB_min : float
        the minimal value of the committor to be accepted as in the ensemble
    pB_max : float
        the maximal value of the committor to be accepted as in the ensemble
    """

    def __init__(
        self,
        trajectories,
        stateA,
        stateB,
        engine,
        storage,
        randomizer,
        n_per_snapshot=20,
        pB_min=0.4,
        pB_max=0.6
    ):
        super(TransitionStateEnsemble, self).__init__(storage)

        # Checks is trajectories is an actual list of trajectories
        try:
            for i in trajectories[0]:
                pass
        except TypeError:
            trajectories = [trajectories]

        self.trajectories = trajectories
        self.stateA = stateA
        self.stateB = stateB
        self.engine = engine
        self.n_per_snapshot = n_per_snapshot
        self.storage = storage
        self.pB_min = pB_min
        self.pB_max = pB_max
        self.pB = 0.0
        self.snap_frame = None
        self.randomizer = randomizer
        self.transition_snapshots = {}
        self.all_steps = []
        # Checks if the pB_min and pB_max input is sane
        if pB_min > pB_max:
            raise ValueError("pB_min should be smaller than pB_max.\n" +
                             "pB_min: " + str(pB_min) + "\n" +
                             "pB_max: " + str(pB_max))

    def next_frame(self, pB=0.0, pB_min=0, pB_max=1,
                   snap_min=None, snap_frame=None, snap_max=None):
        """
        Calculates the next shooting frame index from the committor using the
        bisection method.

        Parameters
        ----------
        pB : float
            The committor of the current snap_frame
        pB_min : float
            The minimal committor of snap_frame to be accepted as
            in the transition state ensemble
        pB_max : float
            The maximum committor of snap_frame to be accepted as
            in the transition state ensemble
        snap_min : int
            The new value of snap_frame must be greater than this value
        snap_frame : int
            The index of the snapshot for which the committor is given
        snap_max : int
            The new value of snap_frame must be smaller than this value
        """

        # Try to get the variables when called without
        if snap_min is None:
            snap_min = self.snap_min
        if snap_frame is None:
            snap_frame = self.snap_frame
        if snap_max is None:
            snap_max = self.snap_max

        # Check if pB limits make sense
        if pB_min > pB_max:
            raise ValueError("pB_min should be smaller than or equal to \
                             pB_min")

        self.output_stream.write("Committor is: "+str(self.pB)+"\n")

        # If the pB is within the interval set snap_frame to None
        if pB_min <= pB <= pB_max:
            self.output_stream.write("Accepted\n\n")
            self.output_stream.flush()
            self.snap_frame = None
            return

        # Shift snap_frame with the bisection method
        if pB < pB_min:
            self.output_stream.write("Moving to later frame\n")
            snap_min = snap_frame
        if pB > pB_max:
            self.output_stream.write("Moving to earlier frame\n")
            snap_max = snap_frame

        self.output_stream.flush()
        snap_frame = snap_min+((snap_max-snap_min)/2)

        # Stops when snap_frame has an invalid value
        # Equal to min or max stop to prevent invinite loops
        if snap_frame <= snap_min or snap_frame >= snap_max:
            if snap_min > snap_max:
                raise ValueError("snap_min should be smaller than snap_max")
            else:
                raise RuntimeError("Shooting frame is outside of shooting " +
                                   "range, you might want to increase " +
                                   "'n_per_snapshot'.\n")

        # Sets the new values
        else:
            self.snap_frame = snap_frame
            self.snap_min = snap_min
            self.snap_max = snap_max

    def run(self):
        """
        Runs the transition state ensemble simulation and returns
        {snap : committor} for each snapshot that is in the ensemble.

        Returns
        -------
        dict :
            This has the shape of {snapshot : committor}
        """

        # Set up variables that are independent of the trajecotory
        randomizer = self.randomizer
        stateA = self.stateA
        stateB = self.stateB

        # Loop over all input trajectories
        for trajectory in self.trajectories:
            self.snap_min = 0
            self.snap_max = len(trajectory)-1

            # Resets snap_frame after 1 loop
            if self.snap_frame is None:
                self.snap_frame = ((self.snap_max-self.snap_min)/2)

            # Loops while there is a next selected frame
            while self.snap_frame is not None:
                self.output_stream.write("Shooting from frame:" +
                                         str(self.snap_frame)+"\n")
                self.output_stream.flush()
                snap = trajectory[self.snap_frame]

                # Sets up the committor simulation
                simulation = paths.CommittorSimulation(storage=self.storage,
                                                       engine=self.engine,
                                                       states=[stateA, stateB],
                                                       randomizer=randomizer,
                                                       direction=None,
                                                       initial_snapshots=[snap]
                                                       )
                # Silences the committor simulation
                devnull = open(os.devnull, 'w')
                simulation.allow_refresh = False
                simulation.output_stream = devnull

                # Runs the committor simulation
                simulation.run(n_per_snapshot=self.n_per_snapshot)
                devnull.close()

                # Append the newly generated steps to self.all_steps

                self.all_steps += self.storage.steps[-self.n_per_snapshot:]
                # Get the shootingpoint analysis object
                # and catch the AssertionError if a trajectory doesn't
                # end in a state
                try:
                    spa = paths.ShootingPointAnalysis(
                              steps=self.storage.steps[-self.n_per_snapshot:],
                              states=[stateA, stateB])
                except AssertionError:
                    raise RuntimeError("Not all of the tries entered a state," +
                                       "you might need to increase" +
                                       "n_frames_max of the engine.")

                # Get the committor from the ShootingPointAnalysis object
                key = spa.step_key(self.storage.steps[-self.n_per_snapshot])

                self.pB = spa.committor(stateB)[key]

                # Get the next shooting frame
                self.next_frame(pB=self.pB,
                                pB_min=self.pB_min,
                                pB_max=self.pB_max)

            # If the snapshot is accepted no next frame is given and the
            # snapshot is stored with its committor value
            self.transition_snapshots[snap] = self.pB

        return self.transition_snapshots
