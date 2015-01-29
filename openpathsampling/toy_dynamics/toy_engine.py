import numpy as np
from openpathsampling.snapshot import Snapshot, Momentum, Configuration
from openpathsampling.dynamics_engine import DynamicsEngine

def convert_to_3Ndim(v):
    ndofs = len(v)
    n_whole_atoms = ndofs / 3
    nrest = 3 - (ndofs % 3)

    out = []
    for i in range(n_whole_atoms):
        out.append([v[3*i+0], v[3*i+1], v[3*i+2]])

    last=[]
    for i in range(ndofs % 3):
        last.append(v[3*n_whole_atoms+i])
    last += [0.0]*nrest

    out.append(last)
    return np.array(out)


def count_atoms(ndofs):
    # first part gives whole atoms, second part says if a partial exists
    return (ndofs / 3) + min(1, ndofs % 3)

class ToyEngine(DynamicsEngine):
    '''The trick is that we have various "simulation" classes (either
    generated directly as here, or subclassed for more complication
    simulation objects as in OpenMM), but they all quack the same when it
    comes to things the DynamicsEngine calls on them for'''

    default_options = { 'pes' : None,
                      'integ' : None,
                      'ndim' : 2,
                      'n_frames_max' : 5000,
                      'nsteps_per_frame' : 10
    }

    def __init__(self, options=None, mode='auto'):
        if 'ndim' not in options:
            options['ndim'] = 2
        options['n_atoms'] = count_atoms(options['ndim'])

        super(ToyEngine, self).__init__(
                                        options=options)

    @property
    def nsteps_per_frame(self):
        return self._nsteps_per_frame

    @nsteps_per_frame.setter
    def nsteps_per_frame(self, value):
        self._nsteps_per_frame = value

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, value):
        self._mass = value
        self.minv = np.reciprocal(value)

    @property
    def current_snapshot(self):
        snap_pos = convert_to_3Ndim(self.positions)
        snap_vel = convert_to_3Ndim(self.velocities)
        snap_pot = self.pes.V(self)
        snap_kin = self.pes.kinetic_energy(self)
        return Snapshot(coordinates=snap_pos,
                        potential_energy=snap_pot,
                        box_vectors=None,
                        velocities=snap_vel,
                        kinetic_energy=snap_kin
                       )

    @current_snapshot.setter
    def current_snapshot(self, snap):
        coords = np.copy(snap.coordinates)
        vels = np.copy(snap.velocities)
        self.positions = np.ravel(coords)[:self.ndim]
        self.velocities = np.ravel(vels)[:self.ndim]

    def generate_next_frame(self):
        self.integ.step(self, self.nsteps_per_frame)
        return self.current_snapshot


    # momentum and configuration properties; these may be removed at some
    # point

    @property
    def momentum(self):
        return Momentum(velocities=convert_to_3Ndim(self.velocities),
                        kinetic_energy=self.pes.kinetic_energy(self)
                       )

    @momentum.setter
    def momentum(self, momentum):
        self.velocities = np.ravel(momentum.velocities)[:self.ndim]

    @property
    def configuration(self):
        return Configuration(coordinates=convert_to_3Ndim(self.positions),
                             box_vectors=None,
                             potential_energy=self.pes.V(self)
                            )

    @configuration.setter
    def configuration(self, configuration):
        self.positions = np.ravel(configuration.coordinates)[:self.ndim]
