"""
Stubs and other tricks used across many tests to get make things quack like
a duck.

@author David W.H. Swenson
"""

import os
from functools import wraps

from pkg_resources import resource_filename
from nose.tools import assert_items_equal, assert_equal, assert_in
import numpy as np
import numpy.testing as npt
import simtk.unit as u

from openpathsampling.trajectory import Trajectory
from openpathsampling.snapshot import Snapshot
from openpathsampling.dynamics_engine import DynamicsEngine
from openpathsampling.topology import Topology
import openpathsampling as paths


def make_1d_traj(coordinates, velocities=None, topology=None):
    if velocities is None:
        velocities = [0.0]*len(coordinates)
    traj = []
    for (pos, vel) in zip(coordinates, velocities):
        snap = Snapshot(coordinates=np.array([[pos, 0, 0]]),
                        velocities=np.array([[vel, 0, 0]]),
                        topology=topology
                        )
        traj.append(snap)
    return Trajectory(traj)

def items_equal(truth, beauty):
    assert_equal(len(truth), len(beauty))
    for (t, b) in zip(truth, beauty):
        if t != b:
            return False
    return True

def assert_items_almost_equal(truth, beauty, tol=10e-7):
    for (t,b) in zip(truth, beauty):
        assert_equal( (t-b)<tol, True)


def assert_equal_array_array(truth, beauty):
    for (t_atom, b_atom) in zip(truth, beauty):
        assert_items_equal(t_atom, b_atom)

def assert_not_equal_array_array(list_a, list_b):
    exist_diff = False
    for (alpha, beta) in zip(list_a, list_b):
        for (elem_a, elem_b) in zip(alpha, beta):
            if elem_a != elem_b:
                exist_diff = True
    return exist_diff

def assert_same_items(list_a, list_b):
    assert_equal(len(list_a), len(list_b))
    for elem_a in list_a:
        assert_in(elem_a, list_b)


class MoverWithSignature(paths.PathMover):
    def __init__(self, input_ensembles, output_ensembles):
        self._in_ensembles = input_ensembles
        self._out_ensembles = output_ensembles

    def move(self, globalstate):
        # need to implement a fake move or this class will be considered abstract
        pass

class CalvinistDynamics(DynamicsEngine):
    def __init__(self, predestination):
        topology = Topology(n_atoms=1, n_spatial=1)
        template = Snapshot(topology=topology)

        super(CalvinistDynamics, self).__init__(options={'n_frames_max' : 12},
                                                template=template)
        self.predestination = make_1d_traj(coordinates=predestination,
                                           velocities=[1.0]*len(predestination),
                                           topology=topology
                                          )
        self.frame_index = None

    @property
    def current_snapshot(self):
        return self._current_snap

    @current_snapshot.setter
    def current_snapshot(self, snap):
        self._current_snap = snap.copy()

    def generate_next_frame(self):
        # find the frame in self.predestination that matches this frame
        if self.frame_index is None:
            for frame in self.predestination:
                frame_val = frame.coordinates[0][0]
                snap_val = self._current_snap.coordinates[0][0]
                # print "looking for " + str(snap_val) + " (" + str(frame_val) + ") " + str(snap_val==frame_val)
                if frame_val == snap_val:
                    self.frame_index = self.predestination.index(frame)
                    break

            #print self.frame_index

        if self._current_snap.velocities[0][0] >= 0:
            self._current_snap = self.predestination[self.frame_index+1].copy()
            self.frame_index += 1
        else:
            self._current_snap = self.predestination[self.frame_index-1].reversed
            self.frame_index -= 1

        return self._current_snap

    def stop(self, trajectory):
        self.frame_index = None

class CallIdentity(object):
    '''Stub for a callable that returns itself'''
    def __init__(self):
        self.name = "Id"
    def __call__(self, value):
        return value


class AtomCounter(object):
    '''Let's be honest: that's all we're using the simulation.system object
    for. So I'll duck-punch.'''
    def __init__(self, n_atoms):
        self.n_atoms = n_atoms

    def getNumParticles(self):
        '''QUAAAAACK'''
        return self.n_atoms

class SimulationDuckPunch(object):
    '''This is what happens when you find a stranger in the Alps.'''
    def __init__(self, topology, system):
        self.system = system
        self.topology = topology

def prepend_exception_message(e, failmsg):
    """Modifies an exception by prepending failmsg"""
    if not e.args:
        e.args = [failmsg]
    else:
        arg0 = failmsg+e.args[0]
        e.args = tuple([arg0] + list(e.args[1:]))

def data_filename(fname, subdir='test_data'):
    return resource_filename('openpathsampling',
                             os.path.join('tests', subdir, fname))

def true_func(value, *args, **kwargs):
    return True

def setify_ensemble_signature(sig):
    return (set(sig[0]), set(sig[1]))


def reorder_ensemble_signature(sig, match_with):
    setified = setify_ensemble_signature(sig)
    found_sigs = []
    for s in match_with:
        if setified == setify_ensemble_signature(s):
            found_sigs.append(s)
    if len(found_sigs) == 0:
        raise RuntimeError("Signature not found for matching: " + repr(sig))
    elif len(found_sigs) > 1:
        raise RuntimeError("More than one form found for signature: " +
                           repr(sig) + "\n" + repr(found_sigs))
    else:
        return found_sigs[0]

def assert_close_unit(v1, v2, *args, **kwargs):
    if type(v1) is u.Quantity:
        assert(v1.unit == v2.unit)
        npt.assert_allclose(v1._value, v2._value, *args, **kwargs)
    else:
        npt.assert_allclose(v1, v2, *args, **kwargs)

def compare_snapshot(snapshot1, snapshot2):
    assert_close_unit(snapshot1.box_vectors, snapshot2.box_vectors, rtol=1e-7, atol=0)
    assert_close_unit(snapshot1.coordinates, snapshot2.coordinates, rtol=1e-7, atol=0)
    assert_close_unit(snapshot1.velocities, snapshot2.velocities, rtol=1e-7, atol=0)

    assert_equal(snapshot1.potential_energy, snapshot2.potential_energy)
    assert_equal(snapshot1.kinetic_energy, snapshot2.kinetic_energy)

class RandomMDEngine(paths.DynamicsEngine):
    _default_options = {}

    def __init__(self, template=None):
        self.options = {
        }

        super(RandomMDEngine, self).__init__(
            options={},
            template=template
        )

        self.initialized = True

    def _build_current_snapshot(self):
        # TODO: Add caching for this and mark if changed

        tmp = self.template

        coordinates = u.Quantity(
            tmp.coordinates._value + np.random.normal(0.0, 0.02, tmp.coordinates.shape),
            tmp.coordinates.unit)
        velocities = u.Quantity(
            np.random.normal(0.0, 0.02, tmp.velocities.shape),
            tmp.velocities.unit)

        return paths.Snapshot(coordinates = coordinates,
                        box_vectors = tmp.box_vectors,
                        potential_energy = tmp.potential_energy,
                        velocities = velocities,
                        kinetic_energy = tmp.kinetic_energy,
                        topology = self.topology
                       )

    @property
    def current_snapshot(self):
        if self._current_snapshot is None:
            self._current_snapshot = self._build_current_snapshot()

        return self._current_snapshot

    @current_snapshot.setter
    def current_snapshot(self, snapshot):
        self._current_snapshot = snapshot

    def generate_next_frame(self):
        self._current_snapshot = None
        return self.current_snapshot

def raises_with_message_like(err, message=None):
    """
    Decorator that allows to run nosetests with raises and testing if the message starts with a txt.

    Notes
    -----
    We use this to check for abstract classes using
    >>> @raises_with_message_like(TypeError, "Can't instantiate abstract class")
    """
    def decorator(fnc):

        @wraps(fnc)
        def _wrapper(*args, **kwargs):
            try:
                fnc(*args, **kwargs)
            except err as e:
                if message is not None and not str(e).startswith(message):
                    raise

        return _wrapper

    return decorator
