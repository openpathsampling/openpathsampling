import numpy as np
import random
import openpathsampling as paths
import openpathsampling.engines.toy as toy
import openpathsampling.high_level.move_strategy as strategies
from openpathsampling import EngineMover
from openpathsampling.high_level.move_strategy import MoveStrategy
from openpathsampling.high_level.move_strategy import levels
from openpathsampling.shooting import NullShootingSelector


class ShiftingMover(EngineMover):
    """
    ShiftingMover propagates a path backward and forward by a fixed length
    and cuts the trajectory in the opposite direction of propagation. The
    propagation direction is chosen from a random uniform distribution.

    Parameters
    ----------
    ensemble : :class:`openpathsampling.Ensemble`
        Ensemble for this shifting mover
    engine ::class:`openpathsampling.DynamicsEngine`
        Engine for the shifting mover
    shift_length ::int
        Length of propagation of the trajectory (default is 3)
    """

    def __init__(self, ensemble, shift_length=3, engine=None):
    #shift length is a time displacement delT from a symmetric distribution around 0
        super(ShiftingMover, self).__init__(
            ensemble=ensemble,
            target_ensemble=ensemble,
            selector=NullShootingSelector(), #fix
	    engine=engine
        )

        self.shift_length = shift_length

    @property
    def direction(self):
        return 'bidirectional' # pragma: no cover


    def _run(self, trajectory , shooting_index):

        shift = random.choice([-1,1]*self.shift_length)
        details = {'Shift': self.shift}

        if shift > 0: #forward shift condition
            #print self.shift_length
            #print 'Initial trajectory', trajectory[1:]

            chopped_traj = trajectory[self.shift_length:] #chop out first shift_length frames

            #generate new trial part of path starting from the end of the chopped traj for fwd shift
            self.engine.current_snapshot = chopped_traj[-1]  # pick last frame
            right_part = self.engine.generate_n_frames(self.shift_length) #gen time steps in future to add to
                                                                          #chopped traj

            #print 'Right part', right_part[1:]
            return  chopped_traj + right_part, details


        elif shift < 0: #backward shift condition
            chopped_traj = trajectory[:-1*self.shift_length]

            #generate new trial part of path starting from the start of the chopped traj for bkwd shift

            self.engine.current_snapshot = chopped_traj[0].reversed  # pick first frame
	    #print self.shift_length, self.engine.current_snapshot

            left_part = self.engine.generate_n_frames(self.shift_length)


            return  left_part.reversed + chopped_traj, details

        #else:
            #print("shift length must be greater or less than 0")


class ShiftingStrategy(MoveStrategy):
    """
    Takes a given scheme and makes the shifting mover. Random choice of backwards or forwards shifting.


    Parameters
    ----------
    shift_length : int
            number of frames in the trajectory to be shifted; default 3
    ensembles : list of :class: `.Ensemble`
        ensembles for which this strategy applies; default None
    engine : :class: `.DynamicsEngine`
        engine to generate the shifted trajectory
    group : str
        The group this strategy is associated with; default "shifting"
    replace : bool
        whether to replace existing movers in the group; default True

    """
    _level = levels.MOVER
    def __init__(self, shift_length=3, ensembles=None, engine=None,
                 group="shifting", replace=True):
        super(ShiftingStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )

        self.shift_length = shift_length
        self.engine = engine

    def make_movers(self, scheme):
        ensemble_list = self.get_ensembles(scheme, self.ensembles)
        ensembles = reduce(list.__add__, map(lambda x: list(x), ensemble_list))
        mover_list=[] #create empty list
        for ensemble in ensembles:
            mover = ShiftingMover(ensemble, self.shift_length, self.engine) #apply the ShiftingMover to each ensemble
            mover_list.append(mover) #add each shifted ensemble to a list

        return mover_list #return out this list




