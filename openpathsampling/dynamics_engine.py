'''
Created on 01.07.2014

@author JDC Chodera
@author: JH Prinz
'''

import logging

import simtk.unit as u

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject

logger = logging.getLogger(__name__)

#=============================================================================
# SOURCE CONTROL
#=============================================================================

__version__ = "$Id: NoName.py 1 2014-07-06 07:47:29Z jprinz $"

#=============================================================================
# Multi-State Transition Interface Sampling
#=============================================================================


class DynamicsEngine(StorableNamedObject):
    '''
    Wraps simulation tool (parameters, storage, etc.)

    Notes
    -----
    Should be considered an abstract class: only its subclasses can be
    instantiated.
    '''

    FORWARD = 1
    BACKWARD = -1

    _default_options = {
        'n_frames_max' : None,
        'timestep' : None
    }

    units = {
        'length' : u.Unit({}),
        'velocity' : u.Unit({}),
        'energy' : u.Unit({})
    }

    def __init__(self, options=None, template=None):
        '''
        Create an empty DynamicsEngine object
        
        Notes
        -----
        The purpose of an engine is to create trajectories and keep track
        of the results. The main method is 'generate' to create a
        trajectory, which is a list of snapshots and then can store the in
        the associated storage. In the initialization this storage is
        created as well as the related Trajectory and Snapshot classes are
        initialized.
        '''

        super(DynamicsEngine, self).__init__()

        self.template = template

        # Trajectories need to know the engine as a hack to get the topology.
        # Better would be a link to the topology directly. This is needed to create
        # mdtraj.Trajectory() objects

        # TODO: Remove this and put the logic outside of the engine. The engine in trajectory is only
        # used to get the solute indices which should depend on the topology anyway
        # Trajectory.engine = self

        self._check_options(options)

        # as default set a newly generated engine as the default engine
        self.set_as_default()

    def _check_options(self, options = None):
        """
        This will register all variables in the options dict as a member variable if
        they are present in either the DynamicsEngine.default_options or this
        classes default_options, no multiple inheritance is supported!
        It will use values with the priority in the following order
        - DynamicsEngine.default_options
        - self.default_options
        - self.options (usually not used)
        - options (function parameter)
        Parameters are only registered if
        1. the variable name is present in the defaults
        2. the type matches the one in the defaults
        3. for variables with units also the units need to be compatible

        Parameters
        ----------
        options : dict of { str : value }
            A dictionary

        Notes
        -----
        Options are what is necessary to recreate the engine, but not runtime variables or independent
        variables like the actual initialization status, the runners or an attached storage.
        If there are non-default options present they will be ignored (no error thrown)
        """
        # start with default options from a dynamics engine
        my_options = {}
        okay_options = {}

        # self.default_options overrides default ones from DynamicsEngine
        for variable, value in self.default_options.iteritems():
            my_options[variable] = value

        if hasattr(self, 'options') and self.options is not None:
            # self.options overrides default ones
            for variable, value in self.options.iteritems():
                my_options[variable] = value

        if options is not None:
            # given options override even default and already stored ones
            for variable, value in options.iteritems():
                my_options[variable] = value

        if my_options is not None:
            for variable, default_value in self.default_options.iteritems():
                # create an empty member variable if not yet present
                if not hasattr(self, variable):
                    okay_options[variable] = None

                if variable in my_options:
                    if type(my_options[variable]) is type(default_value):
                        if type(my_options[variable]) is u.Unit:
                            if my_options[variable].unit.is_compatible(default_value):
                                okay_options[variable] = my_options[variable]
                            else:
                                raise ValueError('Unit of option "' + str(variable) + '" (' + str(my_options[variable].unit) + ') not compatible to "' + str(default_value.unit) + '"')

                        elif type(my_options[variable]) is list:
                            if type(my_options[variable][0]) is type(default_value[0]):
                                okay_options[variable] = my_options[variable]
                            else:
                                raise ValueError('List elements for option "' + str(variable) + '" must be of type "' + str(type(default_value[0])) + '"')
                        else:
                            okay_options[variable] = my_options[variable]
                    elif isinstance(type(my_options[variable]), type(default_value)):
                        okay_options[variable] = my_options[variable]
                    elif default_value is None:
                        okay_options[variable] = my_options[variable]
                    else:
                        raise ValueError('Type of option "' + str(variable) + '" (' + str(type(my_options[variable])) + ') is not "' + str(type(default_value)) + '"')

            self.options = okay_options
        else:
            self.options = {}

    def __getattr__(self, item):
        # default is to look for an option and return it's value
        return self.options[item]

    @property
    def topology(self):
        return self.template.topology

    @property
    def n_atoms(self):
        return self.topology.n_atoms

    @property
    def n_spatial(self):
        return self.topology.n_spatial

    def to_dict(self):
        return {
            'options' : self.options,
            'template' : self.template
        }

    def set_as_default(self):
        paths.EngineMover.engine = self

    @property
    def default_options(self):
        default_options = {}
        default_options.update(DynamicsEngine._default_options)
        default_options.update(self._default_options)
        return default_options

    def start(self, snapshot=None):
        if snapshot is not None:
            self.current_snapshot = snapshot

    def stop(self, trajectory):
        """Nothing special needs to be done for direct-control simulations
        when you hit a stop condition."""
        pass

    def stop_conditions(self, trajectory, continue_conditions=None, 
                        trusted=True):
        """
        Test whether we can continue; called by generate a couple of times,
        so the logic is separated here.

        Parameters
        ----------
        trajectory : Trajectory
            the trajectory we've generated so far
        continue_conditions : list of function(Trajectory)
            callable function of a 'Trajectory' that returns True or False.
            If one of these returns False the simulation is stopped.

        Returns
        -------
        boolean:
            true if the dynamics should be stopped; false otherwise
        """
        stop = False
        if continue_conditions is not None:
            for condition in continue_conditions:
                can_continue = condition(trajectory, trusted)
                stop = stop or not can_continue
        return stop


    def generate_forward(self, snapshot, ensemble):
        """
        Generate a potential trajectory in ensemble simulating forward in time
        """

        return self.generate(snapshot, ensemble.can_append, direction=+1)

    def generate_backward(self, snapshot, ensemble):
        """
        Generate a potential trajectory in ensemble simulating forward in time
        """

        return self.generate(snapshot, ensemble.can_prepend, direction=-1)

    def generate(self, snapshot, running=None, direction=+1):
        r"""
        Generate a trajectory consisting of ntau segments of tau_steps in
        between storage of Snapshots.

        Parameters
        ----------
        snapshot : Snapshot 
            initial coordinates and velocities in form of a Snapshot object
        running : (list of) function(Trajectory)
            callable function of a 'Trajectory' that returns True or False.
            If one of these returns False the simulation is stopped.
        direction : -1 or +1 (DynamicsEngine.FORWARD or DynamicsEngine.BACKWARD)
            If +1 then this will integrate forward, if -1 it will reversed the
            momenta of the given snapshot and then prepending generated snapshots
            with reversed momenta. This will generate a _reversed_ trajectory that
            effectively ends in the initial snapshot

        Returns
        -------    
        trajectory : Trajectory
            generated trajectory of initial conditions, including initial
            coordinate set

        Notes
        -----
        If the returned trajectory has length n_frames_max it can still happen
        that it stopped because of the stopping criterion. You need to check
        in that case.
        """

        if direction == 0:
            raise RuntimeError('direction must be positive (FORWARD) or negative (BACKWARD).')

        try:
            iter(running)
        except:
            running = [running]

        trajectory = paths.Trajectory()

        if direction > 0:
            self.current_snapshot = snapshot
        elif direction < 0:
            # backward simulation needs reversed snapshots
            self.current_snapshot = snapshot.reversed

        self.start()

        # Store initial state for each trajectory segment in trajectory.
        trajectory.append(snapshot)

        frame = 0
        # maybe we should stop before we even begin?
        stop = self.stop_conditions(trajectory=trajectory,
                                    continue_conditions=running,
                                    trusted=False)

        logger.info("Starting trajectory")
        log_freq = 10 # TODO: set this from a singleton class
        while stop == False:
            if self.options.get('n_frames_max', None) is not None :
                if len(trajectory) >= self.options['n_frames_max']:
                    break

            # Do integrator x steps
            snapshot = self.generate_next_frame()
            frame += 1
            if frame % log_freq == 0:
                logger.info("Through frame: %d", frame)

            # Store snapshot and add it to the trajectory. Stores also
            # final frame the last time
            if direction > 0:
                trajectory.append(snapshot)
            elif direction < 0:
                # We are simulating forward and just build in backwards order
                trajectory.prepend(snapshot.reversed)

            # Check if we should stop. If not, continue simulation
            stop = self.stop_conditions(trajectory=trajectory,
                                        continue_conditions=running)

        # exit the while loop once we must stop, so we call the engine's
        # stop function (which should manage any end-of-trajectory
        # cleanup)
        self.stop(trajectory)
        logger.info("Finished trajectory, length: %d", frame)
        return trajectory

    def generate_next_frame(self):
        raise NotImplementedError('Next frame generation must be implemented!')
