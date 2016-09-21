# coding: utf-8# Alanine Multistate

## Example Simulation





### Import and general setup

# standard packages
import numpy as np
import mdtraj as md
import pandas as pd
import math
import random
import time

# helpers for phi-psi plotting
import alatools as ala
import matplotlib.pyplot as plt

# OpenMM
from simtk.openmm import app
import simtk.openmm as mm
import simtk.unit as unit

# OpenPathSampling
import openpathsampling as paths
import openmmtools as omt

# Visualization of PathTrees
from openpathsampling.visualize import PathTree
from IPython.display import SVG

# the openpathsampling OpenMM engine
import openpathsampling.engines.openmm as eng
### Set simulation options and create a simulator object



template = eng.snapshot_from_pdb("../data/Alanine_solvated.pdb")
##### 1. the force field

forcefield = app.ForceField('amber96.xml', 'tip3p.xml')
##### 2. the system object

pdb = app.PDBFile("../data/Alanine_solvated.pdb")

system = forcefield.createSystem(
    pdb.topology, 
    nonbondedMethod=app.PME, 
    nonbondedCutoff=1.0*unit.nanometers,
    constraints=app.HBonds, 
    rigidWater=True,
    ewaldErrorTolerance=0.0005
)
##### 3. the integrator

integrator = omt.integrators.VVVRIntegrator(
    temperature=300 * unit.kelvin,         # temperature
    timestep=2.0 * unit.femtoseconds    # integration step size
)
integrator.setConstraintTolerance(0.00001)
##### 4. the platform

##### 5. OpenMM properties

openmm_properties = {'OpenCLPrecision': 'mixed'}
##### 6. OPS options

engine_options = {
    'n_frames_max': 5000,
    'nsteps_per_frame': 10
}
engine = eng.Engine(
    template.topology, 
    system, 
    integrator, 
    openmm_properties=openmm_properties,
    options=engine_options
)
engine.name = 'default'
integrator_high = mm.LangevinIntegrator(
    1000 * unit.kelvin,         # temperature
    1.0 / unit.picoseconds,    # friction coefficient
    1.0 * unit.femtoseconds    # integration step size
)
integrator.setConstraintTolerance(0.00001)
engine_high_options = {
    'n_frames_max': 10000,
    'nsteps_per_frame': 20
}
engine_high = engine.from_new_options(
    integrator=integrator_high,
    options=engine_high_options)
engine_high.name = 'high'
paths.EngineMover.engine = engine
engine.initialize()
engine_high.initialize()
### Equilibrate

# engine_high.current_snapshot = template
# engine_high.minimize()
# initial_snapshot_high = engine_high.current_snapshot
### Create the storage

storage = paths.Storage("ala_mstis_bootstrap.nc", 'w')
storage.save(engine);
storage.save(engine_high);
storage.tag['template'] = template
### State Definitions

states = ['A', 'B', 'C', 'D', 'E', 'F']
# states = ['A', 'B', 'C', 'D']
state_centers = {
    'A' : [-150, 150], 
    'B' : [-70, 135], 
    'C' : [-150, -65], 
    'D' : [-70, -50], 
    'E' : [50, -100], 
    'F' : [40, 65]
}
interface_levels = {
    'A' : [10, 20, 45, 65, 80],
    'B' : [10, 20, 45, 65, 75],
    'C' : [10, 20, 45, 60],
    'D' : [10, 20, 45, 60],
    'E' : [10, 20, 45, 65, 80],
    'F' : [10, 20, 45, 65, 80],
}
storage.tag['states'] = states
storage.tag['state_centers'] = state_centers
storage.tag['interface_levels'] = interface_levels
### Order Parameters





psi_atoms = [6,8,14,16]
psi = paths.MDTrajFunctionCV(
    name="psi", 
    f=md.compute_dihedrals,
    topology=template.topology,
    indices=[psi_atoms]
).with_diskcache()

phi_atoms = [4,6,8,14]
phi = paths.MDTrajFunctionCV(
    name="phi", 
    f=md.compute_dihedrals,
    topology=template.topology,
    indices=[phi_atoms]
).with_diskcache()

storage.save([psi, phi]);
def circle_degree(snapshot, center, phi, psi):
    import math
    degrees = 180/3.14159    
    psi_deg = psi(snapshot) * degrees
    phi_deg = phi(snapshot) * degrees
    return math.sqrt(
        min( phi_deg - center[0], 360 - phi_deg + center[0])**2 + 
        min( psi_deg - center[1], 360 - psi_deg + center[1])**2)
cv_state = dict()
for state in state_centers:
    op = paths.FunctionCV(
        name = 'op' + state,
        f=circle_degree,
        center=state_centers[state],
        psi=psi,
        phi=phi
    )
    cv_state[state] = op
### Volumes

interface_sets = {}
for state, levels in interface_levels.iteritems():
    interface_sets[state] = \
        paths.VolumeInterfaceSet(cv_state[state], 0.0, levels)
vol_state = {}
for state, levels in interface_levels.iteritems():
    vol_state[state] = interface_sets[state][0]
    vol_state[state].name = state
### Visualize in Phi/Psi space

# reload(ala)
plot = ala.TwoCVSpherePlot(
    cvs=(phi, psi),
    states=[vol_state[vol] for vol in states],
    state_centers=[state_centers[vol] for vol in states],
    interface_levels=[interface_levels[vol] for vol in states]
)

plot.new()
# plot the phi/psi plot and show all states and interfaces
plot.main()
# add the initial template
plot.add_snapshot(template, label='initial')
### Set up the MSTIS network

ms_outers = paths.MSOuterTISInterface.from_lambdas(
    {ifaces: max(ifaces.lambdas) + 5
     for state, ifaces in interface_sets.items()}
)
mstis = paths.MSTISNetwork([
    (vol_state[state], interface_sets[state])          # core, interface set
    for state in states],
    ms_outers=ms_outers
)
storage.tag['network'] = mstis
### Set up the `MoveScheme`



scheme = paths.DefaultScheme(mstis)
### Initial trajectories

hit_all_states_emsemble = paths.join_ensembles([paths.AllOutXEnsemble(vol_state[state]) for state in states])
# initial_trajectories = engine_high.generate(template, [hit_all_states])
trajectory = paths.Trajectory([template])
# generate the iterator starting with the current (old) trajectory
it = engine_high.iter_generate(
    trajectory, 
    [hit_all_states_emsemble], 
    intervals=100)

# run the iterator until it is finished
for traj in it:
    # get a list of found states for status report
    found_states = ','.join(
        state for state in states 
        if paths.PartInXEnsemble(vol_state[state])(traj))

    # output status report
    s = 'Ran %6d steps [%s]. Found states [%s]' % (
        len(traj) - 1, 
        (len(traj) - 1) * engine.snapshot_timestep,
        found_states)
    paths.tools.refresh_output(s)
    
    # remember the last output for later. Python 2 will leak `traj` into
    # local, but Python 3 does not and we want to be compatible
    trajectory = traj
plot.new()
plot.main()
plot.add_trajectory(trajectory)
### Find initial samples

























total_sample_set = scheme.initial_conditions_from_trajectories(
    trajectory, 
#     sample_set = total_sample_set,
    # this strategy is useful for MSTIS Networks and similar ones when generating
    # from a single long trajectory that could contain MinusInterface samples
    strategies = [
        # 1. split and pick shortest and exclude for MinusInterfaceEnsembles
        ('split', {'unique': 'shortest', 'exclude': paths.MinusInterfaceEnsemble}),  
        # 2. split and pick median length
        ('split', {'unique': 'median'}),
        # 3. extend-complex (implemented for minus with segments A-X-A)
        'extend-complex',
        # 4. extend-minimal (implemented for minus and tis with crossings A-X)
        'extend-minimal'],
    engine=engine)
print scheme.initial_conditions_report(total_sample_set)
# loop until we are done
while scheme.check_initial_conditions(total_sample_set)[0]:
    total_sample_set = scheme.initial_conditions_from_trajectories(
        trajectory, 
        sample_set = total_sample_set,
        strategies = [
            'extend-complex',
            'extend-minimal'],
        engine=engine)
## Equilibration





equil_scheme = paths.MoveScheme(mstis)

# tell the scheme to actually use OneWayShooting and nothing else
import openpathsampling.analysis.move_strategy as strat
equil_scheme.append([
        strat.OneWayShootingStrategy(engine=engine), 
        strat.OrganizeByMoveGroupStrategy()
    ])
equilibration = paths.PathSampling(
    storage=None,
    sample_set=total_sample_set,
    move_scheme=equil_scheme,
)
equilibration.run(250)
equilibrated_sset = equilibration.sample_set
engine_list = {}

unique_snapshots = set(sum(
    [
        samp.trajectory for samp in equilibrated_sset 
        if isinstance(samp.ensemble, paths.TISEnsemble)
    ], 
    []))

for samp in equilibrated_sset:
    if isinstance(samp.ensemble, paths.TISEnsemble):
        for snap in samp.trajectory:
            eng = snap.engine
            engine_list[eng] = engine_list.get(eng, 0) + 1
        
for eng, counts in engine_list.items():
    print eng.name, counts
storage.tag['sampleset'] = equilibrated_sset.copy_without_parents()
plt.figure(21, (12,5))
plt.subplot(121)
plt.title('High temperature')
plot.main()
for s in total_sample_set:
    plot.add_trajectory(s.trajectory, line=True)
    
plt.subplot(122)
plt.title('Equilibrated')
plot.main()
plot.zoom = 180/3.1415926
for s in equilibrated_sset:
    plot.add_trajectory(s.trajectory, line=True)
equilibrated_sset.sanity_check()
print 'snapshots:', len(storage.snapshots)
print 'trajectories:', len(storage.trajectories)
print 'samples:', len(storage.samples)
print 'filesize:', storage.file_size_str
# storage.close()

