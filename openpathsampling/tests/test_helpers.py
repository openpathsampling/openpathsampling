"""
Stubs and other tricks used across many tests to get make things quack like
a duck.

@author David W.H. Swenson
"""

import os
from pkg_resources import resource_filename 
from nose.tools import assert_items_equal, assert_equal

from openpathsampling.trajectory import Trajectory
from openpathsampling.snapshot import Snapshot
from openpathsampling.dynamics_engine import DynamicsEngine
from openpathsampling.topology import Topology
import numpy as np

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

            print self.frame_index

        if self._current_snap.velocities[0][0] >= 0:
            self._current_snap = self.predestination[self.frame_index+1].copy()
            self.frame_index += 1
        else:
            self._current_snap = self.predestination[self.frame_index-1].reversed_copy()
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
    def __init__(self, natoms):
        self.natoms = natoms

    def getNumParticles(self):
        '''QUAAAAACK'''
        return self.natoms

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

def true_func(value):
    return True
