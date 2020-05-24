"""
Created on 01.07.2014

@author JDC Chodera
@author: JH Prinz
"""

import logging
import sys

from openpathsampling.netcdfplus import StorableNamedObject
from openpathsampling.integration_tools import is_simtk_unit_type

from .snapshot import BaseSnapshot
from .trajectory import Trajectory

from .delayedinterrupt import DelayedInterrupt

logger = logging.getLogger(__name__)

if sys.version_info > (3, ):
    basestring = str

# =============================================================================
# SOURCE CONTROL
# =============================================================================

__version__ = "$Id: NoName.py 1 2014-07-06 07:47:29Z jprinz $"


# =============================================================================
# Base dynamics engine class
# =============================================================================

class EngineError(Exception):
    def __init__(self, message, last_trajectory):
        # Call the base class constructor with the parameters it needs
        super(EngineError, self).__init__(message)

        # Now for your custom code...
        self.last_trajectory = last_trajectory


class EngineMaxLengthError(EngineError):
    pass


class EngineNaNError(EngineError):
    pass


class DynamicsEngine(StorableNamedObject):
    """
    Wraps simulation tool (parameters, storage, etc.)

    Attributes
    ----------
    on_nan : str
        set the behaviour of the engine when `NaN` is detected.
        Possible is

        1.  `fail` will raise an exception `EngineNaNError`
        2.  `retry` will rerun the trajectory in engine.generate, these moves
            do not satisfy detailed balance

    on_error : str
        set the behaviour of the engine when an exception happens.
        Possible is

        1.  `fail` will raise an exception `EngineError`
        2.  `retry` will rerun the trajectory in engine.generate, these moves
            do not satisfy detailed balance

    on_max_length : str
        set the behaviour if the trajectory length is `n_frames_max`.
        If `n_frames_max == 0` this will be ignored and nothing happens.
        Possible is

        1.  `fail` will raise an exception `EngineMaxLengthError`
        2.  `stop` will stop and return the max length trajectory (default)
        3.  `retry` will rerun the trajectory in engine.generate, these moves
            do not satisfy detailed balance

    retries_when_nan : int, default: 2
        the number of retries (if chosen) before an exception is raised

    retries_when_error : int, default: 2
        the number of retries (if chosen) before an exception is raised

    retries_when_max_length : int, default: 0
        the number of retries (if chosen) before an exception is raised

    on_retry : str or callable
        the behaviour when a try is started. Since you have already generated
        some trajectory you might not restart completely. Possibilities are

        1.  `full` will restart completely and use the initial frames (default)
        2.  `keep_half` will cut the existing in half but keeping at least the initial
        3.  `remove_interval` will remove as many frames as the `interval`
        4.  a callable will be used as a function to generate the new from the
            old trajectories, e.g. `lambda t: t[:10]` would restart with the
            first 10 frames

    Notes
    -----
    Should be considered an abstract class: only its subclasses can be
    instantiated.
    """

    FORWARD = 1
    BACKWARD = -1

    _default_options = {
        'n_frames_max': None,
        'on_max_length': 'fail',
        'on_nan': 'fail',
        'retries_when_nan': 2,
        'retries_when_error': 0,
        'retries_when_max_length': 0,
        'on_retry': 'full',
        'on_error': 'fail'
    }

    #units = {
        #'length': u.Unit({}),
        #'velocity': u.Unit({}),
        #'energy': u.Unit({})
    #}
    units = {
        'length': None,
        'velocity': None,
        'energy': None
    }

    base_snapshot_type = BaseSnapshot

    def __init__(self, options=None, descriptor=None):
        """
        Create an empty DynamicsEngine object

        Notes
        -----
        The purpose of an engine is to create trajectories and keep track
        of the results. The main method is 'generate' to create a
        trajectory, which is a list of snapshots and then can store the in
        the associated storage. In the initialization this storage is
        created as well as the related Trajectory and Snapshot classes are
        initialized.
        """

        super(DynamicsEngine, self).__init__()

        self.descriptor = descriptor
        self._check_options(options)
        self.interrupter = DelayedInterrupt

    @property
    def current_snapshot(self):
        return None

    @current_snapshot.setter
    def current_snapshot(self, snap):
        pass

    def to_dict(self):
        return {
            'options': self.options,
            'descriptor': self.descriptor
        }

    def _check_options(self, options=None):
        """
        This will register all variables in the options dict as a member
        variable if they are present in either the
        `DynamicsEngine.default_options` or this
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
        Options are what is necessary to recreate the engine, but not runtime
        variables or independent variables like the actual initialization
        status, the runners or an attached storage.
        If there are non-default options present they will be ignored
        (no error thrown)
        """
        # start with default options from a dynamics engine
        my_options = {}
        okay_options = {}

        # self.default_options overrides default ones from DynamicsEngine
        for variable, value in self.default_options.items():
            my_options[variable] = value

        if hasattr(self, 'options') and self.options is not None:
            # self.options overrides default ones
            for variable, value in self.options.items():
                my_options[variable] = value

        if options is not None:
            # given options override even default and already stored ones
            for variable, value in options.items():
                my_options[variable] = value

        if my_options is not None:
            for variable, default_value in self.default_options.items():
                # create an empty member variable if not yet present
                if not hasattr(self, variable):
                    okay_options[variable] = None

                if variable in my_options:
                    if type(my_options[variable]) is type(default_value):
                        #if type(my_options[variable]) is u.Unit:
                        if is_simtk_unit_type(my_options[variable]):
                            if my_options[variable].unit.is_compatible(
                                    default_value):
                                okay_options[variable] = my_options[variable]
                            else:
                                raise ValueError(
                                    'Unit of option "' + str(variable) + '" (' +
                                    str(my_options[variable].unit) +
                                    ') not compatible to "' +
                                    str(default_value.unit) +
                                    '"')

                        elif type(my_options[variable]) is list:
                            if isinstance(
                                    my_options[variable][0],
                                    type(default_value[0])):
                                okay_options[variable] = my_options[variable]
                            else:
                                raise \
                                    ValueError(
                                        'List elements for option "' +
                                        str(variable) + '" must be of type "' +
                                        str(type(default_value[0])) + '"')
                        else:
                            okay_options[variable] = my_options[variable]
                    elif isinstance(my_options[variable], type(default_value)):
                        okay_options[variable] = my_options[variable]
                    elif isinstance(my_options[variable], basestring) \
                            and isinstance(default_value, basestring):
                        okay_options[variable] = my_options[variable]
                    elif default_value is None:
                        okay_options[variable] = my_options[variable]
                    else:
                        raise ValueError(
                            'Type of option "' + str(variable) + '" (' +
                            str(type(my_options[variable])) + ') is not "' +
                            str(type(default_value)) + '"')

            self.options = okay_options
        else:
            self.options = {}

    def __getattr__(self, item):
        # first, check for errors that might be shadowed in properties
        if item in self.__class__.__dict__:
            # we should have this attribute
            p = self.__class__.__dict__[item]
            if isinstance(p, property):
                # re-run, raise the error inside the property
                try:
                    result = p.fget(self)
                except:
                    raise
                else:
                    # alternately, trust the fixed result with
                    # return result  # miraculously fixed
                    raise AttributeError(
                        "Unknown problem occurred in property"
                        + str(p.fget.__name__) + ": Second attempt returned"
                        + str(result)
                    )
            # for now, items in dict that fail with AttributeError will just
            # give the default message; to change, add something here like:
            # raise AttributeError("Something went wrong with " + str(item))

        # see, if the attribute is actually a dimension
        if self.descriptor is not None:
            if item in self.descriptor.dimensions:
                return self.descriptor.dimensions[item]

        # fallback is to look for an option and return it's value
        try:
            # extra step here seems to avoid recursion problem in Py3
            option_dict = object.__getattribute__(self, 'options')
            return option_dict[item]
        except KeyError:
            # convert KeyError to AttributeError
            default_msg = "'{0}' has no attribute '{1}'"
            raise AttributeError(
                (default_msg + ", nor does its options dictionary").format(
                    self.__class__.__name__,
                    item
                )
            )

    @property
    def dimensions(self):
        if self.descriptor is None:
            return {}
        else:
            return self.descriptor.dimensions

    def set_as_default(self):
        import openpathsampling as p
        p.EngineMover.engine = self

    @property
    def default_options(self):
        default_options = {}
        default_options.update(DynamicsEngine._default_options)
        default_options.update(self._default_options)
        return default_options

    # def strip_units(self, item):
        # """Remove units and set in the standard unit set for this engine.

        # Each engine needs to know how to do its own unit system. The default
        # assumes there is no unit system.

        # Parameters
        # ----------
        # item : object with units
            # the input with units

        # Returns
        # -------
        # float or iterable
            # the result without units, in the engine's specific unit system
        # """
        # return item

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
        trajectory : :class:`openpathsampling.trajectory.Trajectory`
            the trajectory we've generated so far
        continue_conditions : (list of) function(Trajectory)
            callable function of a 'Trajectory' that returns True or False.
            If one of these returns False the simulation is stopped.
        trusted : bool
            If `True` (default) the stopping conditions are evaluated
            as trusted.

        Returns
        -------
        bool
            true if the dynamics should be stopped; false otherwise
        """
        stop = False
        if continue_conditions is not None:
            if isinstance(continue_conditions, list):
                for condition in continue_conditions:
                    can_continue = condition(trajectory, trusted)
                    stop = stop or not can_continue
            else:
                stop = not continue_conditions(trajectory, trusted)

        return stop

    def generate(self, snapshot, running=None, direction=+1):
        r"""
        Generate a trajectory consisting of ntau segments of tau_steps in
        between storage of Snapshots.

        Parameters
        ----------
        snapshot : :class:`.Snapshot`
            initial coordinates and velocities in form of a Snapshot object
        running : (list of) function(:class:`.Trajectory`)
            callable function of a 'Trajectory' that returns True or False.
            If one of these returns False the simulation is stopped.
        direction : -1 or +1 (DynamicsEngine.FORWARD or DynamicsEngine.BACKWARD)
            If +1 then this will integrate forward, if -1 it will reversed the
            momenta of the given snapshot and then prepending generated
            snapshots with reversed momenta. This will generate a _reversed_
            trajectory that effectively ends in the initial snapshot

        Returns
        -------
        trajectory : :class:`.Trajectory`
            generated trajectory of initial conditions, including initial
            coordinate set

        Notes
        -----
        If the returned trajectory has length n_frames_max it can still happen
        that it stopped because of the stopping criterion. You need to check
        in that case.
        """

        trajectory = None
        it = self.iter_generate(
            snapshot,
            running,
            direction,
            intervals=0,
            max_length=self.options['n_frames_max'])

        for trajectory in it:
            pass

        return trajectory

    def iter_generate(self, initial, running=None, direction=+1,
                      intervals=10, max_length=0):
        r"""
        Return a generator that will generate a trajectory, returning the
        current trajectory in given intervals

        Parameters
        ----------
        initial : :class:`.Snapshot` or :class:`.Trajectory`
            initial coordinates and velocities in form of a Snapshot object
            or a trajectory
        running : (list of) function(:class:`.Trajectory`)
            callable function of a 'Trajectory' that returns True or False.
            If one of these returns False the simulation is stopped.
        direction : -1 or +1 (DynamicsEngine.FORWARD or DynamicsEngine.BACKWARD)
            If +1 then this will integrate forward, if -1 it will reversed the
            momenta of the given snapshot and then prepending generated
            snapshots with reversed momenta. This will generate a _reversed_
            trajectory that effectively ends in the initial snapshot
        intervals : int
            number steps after which the current status is returned. If `0`
            it will run until the end or a keyboard interrupt is detected
        max_length : int
            will limit the simulation length to a number of steps. Default is
            `0` which will run unlimited

        Yields
        ------
        trajectory : :class:`openpathsampling.trajectory.Trajectory`
            generated trajectory of initial conditions, including initial
            coordinate set

        Notes
        -----
        If the returned trajectory has length n_frames_max it can still happen
        that it stopped because of the stopping criterion. You need to check
        in that case.
        """

        if direction == 0:
            raise RuntimeError(
                'direction must be positive (FORWARD) or negative (BACKWARD).')

        try:
            iter(running)
        except TypeError:
            running = [running]

        if hasattr(initial, '__iter__'):
            initial = Trajectory(initial)
        else:
            initial = Trajectory([initial])

        valid = False
        attempt_nan = 0
        attempt_error = 0
        attempt_max_length = 0
        trajectory = initial

        final_error = None
        errors = []

        while not valid and final_error is None:
            if attempt_nan + attempt_error > 1:
                # let's get a new initial trajectory the way the user wants to
                if self.on_retry == 'full':
                    trajectory = initial
                elif self.on_retry == 'remove_interval':
                    trajectory = \
                        trajectory[:max(
                            len(initial),
                            len(trajectory) - intervals)]
                elif self.on_retry == 'keep_half':
                    trajectory = \
                        trajectory[:min(
                            int(len(trajectory) * 0.9),
                            max(
                                len(initial),
                                int(len(trajectory) / 2)))]
                elif hasattr(self.on_retry, '__call__'):
                    trajectory = self.on_retry(trajectory)
            
            """ Case of run dying before first output"""
            if len(trajectory) >= 1:
                if direction > 0:
                    self.current_snapshot = trajectory[-1]
                elif direction < 0:
                    # backward simulation needs reversed snapshots
                    self.current_snapshot = trajectory[0].reversed

            logger.info("Starting trajectory")
            self.start()

            frame = 0
            # maybe we should stop before we even begin?
            stop = self.stop_conditions(trajectory=trajectory,
                                        continue_conditions=running,
                                        trusted=False)

            log_rate = 10
            has_nan = False
            has_error = False

            while not stop:
                if intervals > 0 and frame % intervals == 0:
                    # return the current status
                    logger.info("Through frame: %d", frame)
                    yield trajectory

                elif frame % log_rate == 0:
                    logger.info("Through frame: %d", frame)

                # Do integrator x steps

                snapshot = None

                try:
                    with self.interrupter():
                        snapshot = self.generate_next_frame()

                        # if self.on_nan != 'ignore' and \
                        if not self.is_valid_snapshot(snapshot):
                            has_nan = True
                            break

                except KeyboardInterrupt as e:
                    # make sure we will report the last state for
                    logger.info('Keyboard interrupt. Shutting down simulation')
                    final_error = e
                    break

                except:
                    # any other error we start a retry
                    e = sys.exc_info()
                    errors.append(e)
                    se = str(e).lower()
                    if 'nan' in se and \
                            ('particle' in se or 'coordinates' in se):
                        # this cannot be ignored because we cannot continue!
                        has_nan = True
                        break
                    else:
                        has_error = True
                        break

                frame += 1

                # Store snapshot and add it to the trajectory.
                # Stores also final frame the last time
                if direction > 0:
                    trajectory.append(snapshot)
                elif direction < 0:
                    trajectory.insert(0, snapshot.reversed)

                if 0 < max_length < len(trajectory):
                    # hit the max length criterion
                    on = self.on_max_length
                    del trajectory[-1]

                    if on == 'fail':
                        final_error = EngineMaxLengthError(
                            'Hit maximal length of %d frames.' %
                            self.options['n_frames_max'],
                            trajectory
                        )
                        break
                    elif on == 'stop':
                        logger.info('Trajectory hit max length. Stopping.')
                        # fail gracefully
                        stop = True
                    elif on == 'retry':
                        attempt_max_length += 1
                        if attempt_max_length > self.retries_when_max_length:
                            if self.on_nan == 'fail':
                                final_error = EngineMaxLengthError(
                                    'Failed to generate trajectory without '
                                    'hitting max length after %d attempts' %
                                    attempt_max_length,
                                    trajectory)
                                break

                if stop is False:
                    # Check if we should stop. If not, continue simulation
                    stop = self.stop_conditions(trajectory=trajectory,
                                            continue_conditions=running)

            if has_nan:
                on = self.on_nan
                if on == 'fail':
                    final_error = EngineNaNError(
                        '`nan` in snapshot', trajectory)
                elif on == 'retry':
                    attempt_nan += 1
                    if attempt_nan > self.retries_when_nan:
                        final_error = EngineNaNError(
                            'Failed to generate trajectory without `nan` '
                            'after %d attempts' % attempt_error,
                            trajectory)

            elif has_error:
                on = self.on_nan
                if on == 'fail':
                    final_error = errors[-1][1]
                    del errors[-1]
                elif on == 'retry':
                    attempt_error += 1
                    if attempt_error > self.retries_when_error:
                        final_error = EngineError(
                            'Failed to generate trajectory without `nan` '
                            'after %d attempts' % attempt_error,
                            trajectory)

            elif stop:
                valid = True

            self.stop(trajectory)

        if errors:
            logger.info('Errors occurred during generation :')
            for no, e in enumerate(errors):
                logger.info('[#%d] %s' % (no, repr(e[1])))

        if final_error is not None:
            yield trajectory
            logger.info("Through frame: %d", len(trajectory))
            raise final_error

        logger.info("Finished trajectory, length: %d", len(trajectory))
        yield trajectory

    def generate_next_frame(self):
        raise NotImplementedError('Next frame generation must be implemented!')

    def generate_n_frames(self, n_frames=1):
        """Generates n_frames, from but not including the current snapshot.

        This generates a fixed number of frames at once. If you desire the
        reversed trajectory, you can reverse the returned trajectory.

        Parameters
        ----------
        n_frames : integer
            number of frames to generate

        Returns
        -------
        paths.Trajectory()
            the `n_frames` of the trajectory following (and not including)
            the initial `current_snapshot`
        """
        self.start()
        traj = Trajectory([self.generate_next_frame()
                           for i in range(n_frames)])
        self.stop(traj)
        return traj

    @staticmethod
    def is_valid_snapshot(snapshot):
        """
        Test the snapshot to be valid. Usually not containing nan

        Returns
        -------
        bool : True
            returns `True` if the snapshot is okay to be used
        """
        return True

    @classmethod
    def check_snapshot_type(cls, snapshot):
        if not isinstance(snapshot, cls.base_snapshot_type):
            logger.warning(
                ('This engine is intended for "%s" and derived classes. '
                 'You are using "%s". Make sure that this is intended.') %
                (cls.base_snapshot_type.__name__, snapshot.__class__.__name__)
            )


class NoEngine(DynamicsEngine):
    _default_options = {}

    def __init__(self, descriptor):
        super(NoEngine, self).__init__()
        self.descriptor = descriptor

    def generate_next_frame(self):
        pass
