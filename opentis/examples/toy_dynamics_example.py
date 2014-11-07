'''
@author: David W.H. Swenson

Example code for using the toy_dynamics package. OpenTIS is designed to take
advantage of the tools provided by OpenMM; however, if you want to use it to
interface with some other MD engine, that's also possible! Here we use the
engine in toy_dynamics.

Also, this gives us a nice way to properly check that the stochastic
thermostat is behaving as desired.
'''

from opentis.snapshot import Snapshot
from opentis.ensemble import LengthEnsemble
from opentis.toy_dynamics.toy_simulation import ToySimulation
from opentis.toy_dynamics.toy_integrators import LangevinBAOABIntegrator
import opentis.toy_dynamics.toy_pes as pes
import numpy as np

if __name__ == "__main__":
    nsteps = 1000000
    my_pes = pes.HarmonicOscillator([1.0, 1.0], [1.0, 1.0], [0.0, 0.0])
    my_integ = LangevinBAOABIntegrator(0.002, 0.5, 1.0)
    sim = ToySimulation(my_pes, my_integ)
    sim.nsteps_per_iteration = 50
    sim.mass = np.array([1.0, 1.0])
    snap = Snapshot(coordinates=np.array([0.0, 0.0]),
                    velocities=np.array([0.1, 0.0]))
    sim.init_simulation_with_snapshot(snap)
    x1 = []
    x2 = []
    v1 = []
    v2 = []
    for i in range(nsteps):
        snap = sim.generate_next_frame()
        # sample the information desired to check distributions
        pos = snap.coordinates
        vel = snap.velocities
        x1.append(pos[0])
        x2.append(pos[1])
        v1.append(vel[0])
        v2.append(vel[1])

    nbins = 50
    rrange = (-2.5, 2.5)
    dens = True
    (x1hist, bins) = np.histogram(x1, bins=nbins, range=rrange, density=dens)
    x2hist = np.histogram(x2, bins=nbins, range=rrange, density=dens)[0]
    v1hist = np.histogram(v1, bins=nbins, range=rrange, density=dens)[0]
    v2hist = np.histogram(v2, bins=nbins, range=rrange, density=dens)[0]

    
    for (b, xx1, xx2, vv1, vv2) in zip(bins, x1hist, x2hist, v1hist, v2hist):
        print b, xx1, xx2, vv1, vv2


