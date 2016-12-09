# coding: utf-8

# standard packages
import mdtraj as md

# OpenMM
from simtk.openmm import app
import simtk.unit as unit

# OpenPathSampling
import openpathsampling as paths
import openmmtools as omt

# the openpathsampling OpenMM engine
import openpathsampling.engines.openmm as eng

# file paths
import config as cf

from openpathsampling.tools import refresh_output

# =============================================================================
# ALANINE MULTISTATE
# =============================================================================
print """ALANINE MULTISTATE"""

# Example Simulation


# -----------------------------------------------------------------------------
# Import and general setup
# -----------------------------------------------------------------------------
print """Import and general setup"""


# -----------------------------------------------------------------------------
# Set simulation options and create a simulator object
# -----------------------------------------------------------------------------
print """Set simulation options and create a simulator object"""

template = eng.snapshot_from_pdb(cf.pdb_file)

print """## 1. the force field"""
forcefield = app.ForceField('amber96.xml', 'tip3p.xml')

print """## 2. the system object"""
pdb = app.PDBFile(cf.pdb_file)

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometers,
    constraints=app.HBonds,
    ewaldErrorTolerance=0.0005
)

print """## 3. the integrator"""
integrator_low = omt.integrators.VVVRIntegrator(
    temperature=300 * unit.kelvin,  # temperature
    timestep=2.0 * unit.femtoseconds  # integration step size
)
integrator_low.setConstraintTolerance(0.00001)

print """## 4. the platform"""

print eng.Engine.available_platforms()

print """## 5. OpenMM properties"""

platform = cf.platform

if platform == 'OpenCL':
    openmm_properties = {'OpenCLPrecision': 'mixed'}
elif platform == 'CUDA':
    openmm_properties = {'CUDAPrecision': 'mixed'}
elif platform == 'CPU':
    openmm_properties = {}
else:
    openmm_properties = {}

print """## 6. OPS options"""
engine_low_options = {
    'n_frames_max': 5000,
    'n_steps_per_frame': 10
}
engine_low = eng.Engine(
    template.topology,
    system,
    integrator_low,
    openmm_properties=openmm_properties,
    options=engine_low_options
)
engine_low.name = 'default'

integrator_high = omt.integrators.VVVRIntegrator(
    temperature=1000 * unit.kelvin,  # temperature
    timestep=2.0 * unit.femtoseconds  # integration step size
)
integrator_high.setConstraintTolerance(0.00001)

engine_high_options = {
    'n_frames_max': 2000,
    'n_steps_per_frame': 10  # twice as many steps with half stepsize
}
engine_high = engine_low.from_new_options(
    integrator=integrator_high,
    options=engine_high_options)
engine_high.name = 'high'

engine_high.initialize(platform)

print 'High-Engine uses'
print 'platform `%s`' % engine_high.platform
print 'temperature `%6.2f K' % (
    float(engine_high.integrator.getGlobalVariableByName('kT')) / 0.0083144621)

# -----------------------------------------------------------------------------
# Equilibrate
# -----------------------------------------------------------------------------
print """Equilibrate"""

# engine_high.current_snapshot = template
# engine_high.minimize()
# initial_snapshot_high = engine_high.current_snapshot

# -----------------------------------------------------------------------------
# Create the storage
# -----------------------------------------------------------------------------
print """Create the storage"""

storage = paths.Storage(cf.storage_setup, 'w')
storage.save(engine_low)
storage.save(engine_high)
storage.tag['template'] = template

# -----------------------------------------------------------------------------
# State Definitions
# -----------------------------------------------------------------------------
print """State Definitions"""

states = ['A', 'B', 'C', 'D', 'E', 'F']
# states = ['A', 'B', 'C', 'D']
state_centers = {
    'A': [-150, 150],
    'B': [-70, 135],
    'C': [-150, -65],
    'D': [-70, -50],
    'E': [50, -100],
    'F': [40, 65]
}
interface_levels = {
    'A': [10, 20, 45, 65, 80],
    'B': [10, 20, 45, 65, 75],
    'C': [10, 20, 45, 60],
    'D': [10, 20, 45, 60],
    'E': [10, 20, 45, 65, 80],
    'F': [10, 20, 45, 65, 80],
}
storage.tag['states'] = states
storage.tag['state_centers'] = state_centers
storage.tag['interface_levels'] = interface_levels

# -----------------------------------------------------------------------------
# Order Parameters
# -----------------------------------------------------------------------------
print """Order Parameters"""

psi_atoms = [6, 8, 14, 16]
psi = paths.MDTrajFunctionCV(
    name="psi",
    f=md.compute_dihedrals,
    topology=template.topology,
    indices=[psi_atoms]
).with_diskcache()

phi_atoms = [4, 6, 8, 14]
phi = paths.MDTrajFunctionCV(
    name="phi",
    f=md.compute_dihedrals,
    topology=template.topology,
    indices=[phi_atoms]
).with_diskcache()

storage.save([psi, phi])


def circle_degree(snapshot, center, phi, psi):
    import numpy
    p = numpy.array([psi(snapshot), phi(snapshot)]) / numpy.pi * 180.0
    delta = numpy.abs(center - p)
    delta = numpy.where(delta > 180.0, delta - 360.0, delta)
    return numpy.hypot(delta[0], delta[1])


cv_state = dict()
for state in state_centers:
    op = paths.FunctionCV(
        name='op' + state,
        f=circle_degree,
        center=state_centers[state],
        psi=psi,
        phi=phi
    )
    cv_state[state] = op

# -----------------------------------------------------------------------------
# Volumes
# -----------------------------------------------------------------------------
print """Volumes"""

interface_sets = {}
for state, levels in interface_levels.iteritems():
    interface_sets[state] = \
        paths.VolumeInterfaceSet(cv_state[state], 0.0, levels)
vol_state = {}
for state, levels in interface_levels.iteritems():
    vol_state[state] = interface_sets[state][0]
    vol_state[state].name = state

# -----------------------------------------------------------------------------
# Set up the MSTIS network
# -----------------------------------------------------------------------------
print """Set up the MSTIS network"""

ms_outers = paths.MSOuterTISInterface.from_lambdas(
    {ifaces: max(ifaces.lambdas) + 5
     for state, ifaces in interface_sets.items()}
)
mstis = paths.MSTISNetwork([
                               (vol_state[state], interface_sets[state])
                               # core, interface set
                               for state in states],
                           ms_outers=ms_outers
                           )
storage.tag['network'] = mstis

# -----------------------------------------------------------------------------
# Set up the `MoveScheme`
# -----------------------------------------------------------------------------
print """Set up the `MoveScheme`"""

scheme = paths.DefaultScheme(mstis)

# -----------------------------------------------------------------------------
# Initial trajectories
# -----------------------------------------------------------------------------
print """Initial trajectories"""

hit_all_states_emsemble = paths.join_ensembles(
    [paths.AllOutXEnsemble(vol_state[state]) for state in states])
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

    visited_states = ','.join(
        state for state in states
        if paths.PartInXEnsemble(vol_state[state])(traj[-100:]))

    # output status report
    s = 'Ran %6d steps [%s]. Found states [%s] / visited [%s]\n' % (
        len(traj) - 1,
        (len(traj) - 1) * engine_high.snapshot_timestep,
        found_states,
        visited_states)

    refresh_output(s)

    # remember the last output for later. Python 2 will leak `traj` into
    # local, but Python 3 does not and we want to be compatible
    trajectory = traj

# -----------------------------------------------------------------------------
# Find initial samples
# -----------------------------------------------------------------------------
print """Find initial samples"""

total_sample_set = scheme.initial_conditions_from_trajectories(
    trajectory,
    strategies=[
        # 1. split and pick shortest and exclude for MinusInterfaceEnsembles
        'split',
        # 2. extend-complex (implemented for minus with segments A-X-A)
        'extend-complex',
        # 3. extend-minimal (implemented for minus and tis with crossings A-X)
        'extend-minimal'],
    engine=engine_high)
print scheme.initial_conditions_report(total_sample_set)

# loop until we are done
while scheme.check_initial_conditions(total_sample_set)[0]:
    total_sample_set = scheme.initial_conditions_from_trajectories(
        trajectory,
        sample_set=total_sample_set,
        strategies=[
            'extend-complex',
            'extend-minimal'],
        engine=engine_high)

# -----------------------------------------------------------------------------
# Equilibration
# -----------------------------------------------------------------------------
print """Equilibration"""

engine_high.unload_context()
engine_low.initialize(cf.platform)

refresh_output('Low-Engine uses', refresh=False)
print 'platform `%s`' % engine_low.platform
print 'temperature `%6.2f K' % (
    float(engine_low.integrator.getGlobalVariableByName('kT')) / 0.0083144621)

print

# Equilibration
equil_scheme = paths.OneWayShootingMoveScheme(mstis, engine=engine_low)

equilibration = paths.PathSampling(
    storage=None,
    sample_set=total_sample_set,
    move_scheme=equil_scheme,
)

storage.sync()

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

equilibrated_sset.sanity_check()
print 'snapshots:', len(storage.snapshots)
print 'trajectories:', len(storage.trajectories)
print 'samples:', len(storage.samples)
print 'filesize:', storage.file_size_str

storage.sync()

storage.close()
