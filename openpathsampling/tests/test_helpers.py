"""
Stubs and other tricks used across many tests to get make things quack like
a duck.

@author David W.H. Swenson
"""

import os
from functools import wraps
from openpathsampling.engines import NoEngine
import numpy as np
import numpy.testing as npt

try:
    import simtk.unit as u
except ImportError:
    u = None

from openpathsampling.integration_tools import is_simtk_quantity_type

try:
    import mdtraj as md
except ImportError:
    md = None

from pkg_resources import resource_filename

import openpathsampling as paths
import openpathsampling.engines.openmm as peng
import openpathsampling.engines.toy as toys
from openpathsampling.engines import Topology

from openpathsampling.engines import DynamicsEngine

def make_1d_traj(coordinates, velocities=None, engine=None):
    if velocities is None:
        velocities = [1.0]*len(coordinates)
    try:
        _ = len(velocities)
    except TypeError:
        velocities = [velocities] * len(coordinates)

    if engine is None:
        engine = toys.Engine(
            {},
            toys.Topology(
                n_spatial=3,
                masses=[1.0, 1.0, 1.0], pes=None
            )
        )
    traj = []
    for (pos, vel) in zip(coordinates, velocities):
        snap = toys.Snapshot(
            coordinates=np.array([[pos, 0, 0]]),
            velocities=np.array([[vel, 0, 0]]),
            engine=engine
        )
        traj.append(snap)
    return paths.Trajectory(traj)

def items_equal(truth, beauty):
    assert len(truth) == len(beauty)
    for (t, b) in zip(truth, beauty):
        if t != b:
            return False
    return True

def assert_items_equal(truth, beauty):
    assert len(truth) == len(beauty)
    for (t, b) in zip(truth, beauty):
        assert t == b

def assert_items_almost_equal(truth, beauty, tol=10e-7):
    for (t,b) in zip(truth, beauty):
        assert abs(t-b) - tol < 0.0


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
    assert len(list_a) == len(list_b)
    for elem_a in list_a:
        assert elem_a in list_b


class MoverWithSignature(paths.PathMover):
    def __init__(self, input_ensembles, output_ensembles):
        super(MoverWithSignature, self).__init__()
        self._in_ensembles = input_ensembles
        self._out_ensembles = output_ensembles

    def move(self, sample_set):
        # need to implement a fake move or this class will be considered abstract
        pass

class CalvinistDynamics(DynamicsEngine):
    def __init__(self, predestination):
        self.topology = Topology(n_atoms=1, n_spatial=1)
        # engine = peng.tools.TopologyEngine(topology)

        super(CalvinistDynamics, self).__init__(options={'n_frames_max': 12})
        self.predestination = make_1d_traj(coordinates=predestination,
                                           velocities=[1.0]*len(predestination),
                                           engine=self
                                          )
        self.frame_index = None

    @property
    def current_snapshot(self):
        return self._current_snap

    @current_snapshot.setter
    def current_snapshot(self, snap):
        self._current_snap = snap

    def generate_next_frame(self):
        # find the frame in self.predestination that matches this frame
        if self.frame_index is None:
            for frame in self.predestination:
                frame_val = frame.coordinates[0][0]
                snap_val = self._current_snap.coordinates[0][0]
                # print "looking for " + str(snap_val) + " (" + str(frame_val) + ") " + str(snap_val==frame_val)
                if abs(frame_val - snap_val) < 1e-7:
                    self.frame_index = self.predestination.index(frame)
                    break

            #print self.frame_index

        if self._current_snap.velocities[0][0] >= 0:
            self._current_snap = self.predestination[self.frame_index+1]
            self.frame_index += 1
        else:
            self._current_snap = self.predestination[self.frame_index-1].reversed
            self.frame_index -= 1

        # print self._current_snap.xyz[0][0]
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
    if is_simtk_quantity_type(v1):
        assert(v1.unit == v2.unit)
        npt.assert_allclose(v1._value, v2._value, *args, **kwargs)
    else:
        npt.assert_allclose(v1, v2, *args, **kwargs)

def compare_snapshot(snapshot1, snapshot2, check_reversed=False):
    if hasattr(snapshot1, 'box_vectors') == hasattr(snapshot2, 'box_vectors'):
        if hasattr(snapshot1, 'box_vectors'):
            assert_close_unit(snapshot1.box_vectors, snapshot2.box_vectors, rtol=1e-7, atol=0)
    else:
        raise AttributeError('Snapshots disagree. Only one uses box_vectors')

    assert_close_unit(snapshot1.coordinates, snapshot2.coordinates, rtol=1e-7, atol=0)
    assert_close_unit(snapshot1.velocities, snapshot2.velocities, rtol=1e-7, atol=0)

    if check_reversed:
        compare_snapshot(snapshot1.reversed, snapshot2.reversed, False)
        assert_close_unit(-1.0 * snapshot1.reversed.velocities, snapshot1.velocities, rtol=1e-7, atol=0)
        assert_close_unit(-1.0 * snapshot2.reversed.velocities, snapshot2.velocities, rtol=1e-7, atol=0)
        assert_close_unit(snapshot1.reversed.coordinates, snapshot1.coordinates, rtol=1e-7, atol=0)
        assert_close_unit(snapshot2.reversed.coordinates, snapshot2.coordinates, rtol=1e-7, atol=0)


class RandomMDEngine(DynamicsEngine):
    _default_options = {}

    def __init__(self, template=None):
        self.options = {
        }

        super(RandomMDEngine, self).__init__()

        self.template = template
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

        return peng.Snapshot.construct(coordinates = coordinates,
                        box_vectors = tmp.box_vectors,
                        velocities = velocities,
                        engine=self
                       )

    @property
    def current_snapshot(self):
        if self._current_snapshot is None:
            self._current_snapshot = self._build_current_snapshot()

        return self._current_snapshot

    @current_snapshot.setter
    def current_snapshot(self, snapshot):
        self._current_snapshot = snapshot

    @property
    def snapshot_timestep(self):
        return 1.0

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

def assert_frame_equal(truth, beauty):
    assert len(truth.index) == len(beauty.index)
    assert len(truth.columns) == len(beauty.columns)
    assert set(truth.index) == set(beauty.index)
    assert set(truth.columns) == set(beauty.columns)
    for idx in truth.index:
        for col in truth.columns:
            truth_val = truth.loc[idx, col]
            beauty_val = beauty.loc[idx, col]
            if np.isnan(truth_val):
                assert np.isnan(beauty_val)
            else:
                assert truth_val == beauty_val

def A2BEnsemble(volume_a, volume_b, trusted=True):
    # this is a little replacement for the same name that used to be in
    # EnsembleFactory. It was only used in tests.
    return paths.SequentialEnsemble([
        paths.AllInXEnsemble(volume_a) & paths.LengthEnsemble(1),
        paths.AllOutXEnsemble(volume_a | volume_b),
        paths.AllInXEnsemble(volume_b) & paths.LengthEnsemble(1)
    ])


class NaNEngine(NoEngine):
    def __init__(self, descriptor):
        super(NaNEngine, self).__init__(descriptor=descriptor)

    @staticmethod
    def is_valid_snapshot(snapshot):
        return False
