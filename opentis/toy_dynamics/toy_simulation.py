import numpy as np
from opentis.snapshot import Snapshot
from opentis.Simulator import Simulator

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



class ToySimulation(Simulator):
    '''The trick is that we have various "simulation" classes (either
    generated directly as here, or subclassed for more complication
    simulation objects as in OpenMM), but they all quack the same when it
    comes to things the Simulator calls on them for'''

    def __init__(self, pes, integ, ndim=2):
        self.pes = pes
        self.integ = integ
        self.ndim = ndim

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

    def load_momentum(self, momentum):
        momentum._velocities = convert_to_3Ndim(self.velocities)
        momentum._kinetic_energy = self.pes.kinetic_energy(self)

    def load_configuration(self, configuration):
        configuration._coordinates = convert_to_3Ndim(self.positions)
        configuration._potential_energy = self.pes.V(self)
        configuration._box_vectors = None # toys without PBCs

    def load_snapshot(self, snapshot):
        snapshot.configuration._coordinates = convert_to_3Ndim(self.positions)
        snapshot.configuration._potential_energy = self.pes.V(self)
        snapshot.momentum._velocities = convert_to_3Ndim(self.velocities)
        snapshot.momentum._kinetic_energy = self.pes.kinetic_energy(self)
        snapshot.configuration._box_vectors = None

    @property
    def current_snaphost(self):
        snapshot = Snapshot()
        snapshot.configuration._coordinates = convert_to_3Ndim(self.positions)
        snapshot.configuration._potential_energy = self.pes.V(self)
        snapshot.momentum._velocities = convert_to_3Ndim(self.velocities)
        snapshot.momentum._kinetic_energy = self.pes.kinetic_energy(self)
        snapshot.configuration._box_vectors = None
        return snapshot

    @current_snaphost.setter
    def current_snaphost(self, snap):
        self.positions = np.ravel(snap.configuration.coordinates)[:self.ndim]
        self.velocities = np.ravel(snap.momentum.velocities)[:self.ndim]

    def generate_next_frame(self):
        self.integ.step(self, self.nsteps_per_frame)
        snap = Snapshot()
        self.load_snapshot(snap)
        return snap

    def start(self, snapshot=None):
        if snapshot is not None:
            self.current_snaphost = snapshot


    def stop(self, trajectory):
        pass # pragma: no cover (no need to test this one)

