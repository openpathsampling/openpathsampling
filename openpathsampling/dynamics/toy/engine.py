import numpy as np

from openpathsampling.dynamics import DynamicsEngine
from snapshot import ToySnapshot as Snapshot

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

    base_snapshot_type = Snapshot

    default_options = {
                      'integ' : None,
                      'n_frames_max' : 5000,
                      'nsteps_per_frame' : 10
    }

    def __init__(self, options, template):
        if 'n_spatial' not in options:
            options['n_spatial'] = template.topology.n_spatial

        options['n_atoms'] = 1

        super(ToyEngine, self).__init__(
                                        options=options)

        self.template = template
        self.mass = template.topology.masses
        self._pes = template.topology.pes

        self.current_snapshot = self.template

    @property
    def pes(self):
        return self._pes

    @property
    def nsteps_per_frame(self):
        return self.options['nsteps_per_frame']

    @nsteps_per_frame.setter
    def nsteps_per_frame(self, value):
        self.options['nsteps_per_frame'] = value

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, value):
        self._mass = value
        self._minv = np.reciprocal(value)

    @property
    def snapshot_timestep(self):
        return self.nsteps_per_frame * self.integ.dt

    @property
    def current_snapshot(self):
        snap_pos = self.positions
        snap_vel = self.velocities
        return Snapshot(
            coordinates=np.array([snap_pos]),
            velocities=np.array([snap_vel]),
            topology=self.template.topology
        )

    @current_snapshot.setter
    def current_snapshot(self, snap):
        self.check_snapshot_type(snap)

        coords = np.copy(snap.coordinates)
        vels = np.copy(snap.velocities)
        self.positions = coords[0]
        self.velocities = vels[0]

    def generate_next_frame(self):
        self.integ.step(self, self.nsteps_per_frame)
        return self.current_snapshot
