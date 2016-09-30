# coding: utf-8


# =============================================================================
# ALANINE MULTISTATE
# =============================================================================
print """ALANINE MULTISTATE"""

# Example Simulation


# -----------------------------------------------------------------------------
# Import and general setup
# -----------------------------------------------------------------------------
print """Import and general setup"""

# standard packages
import mdtraj as md

# OpenMM
from simtk.openmm import app
import simtk.openmm as mm
import simtk.unit as unit

# OpenPathSampling
import openpathsampling as paths
import openmmtools as omt

# the openpathsampling OpenMM engine
import openpathsampling.engines.openmm as eng

# file paths
import config as cf

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
    rigidWater=True,
    ewaldErrorTolerance=0.0005
)

print """## 3. the integrator"""
integrator = omt.integrators.VVVRIntegrator(
    temperature=300 * unit.kelvin,  # temperature
    timestep=2.0 * unit.femtoseconds  # integration step size
)
integrator.setConstraintTolerance(0.00001)

print """## 4. the platform"""

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
    1000 * unit.kelvin,  # temperature
    1.0 / unit.picoseconds,  # friction coefficient
    1.0 * unit.femtoseconds  # integration step size
)
integrator.setConstraintTolerance(0.00001)
engine_high_options = {
    'n_frames_max': 2000,
    'nsteps_per_frame': 20
}
engine_high = engine.from_new_options(
    integrator=integrator_high,
    options=engine_high_options)
engine_high.name = 'high'
paths.EngineMover.engine = engine
engine.initialize()
engine_high.initialize()

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
storage.save(engine)
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
    import math
    degrees = 180 / 3.14159
    psi_deg = psi(snapshot) * degrees
    phi_deg = phi(snapshot) * degrees
    return math.sqrt(
        min(phi_deg - center[0], 360 - phi_deg + center[0]) ** 2 +
        min(psi_deg - center[1], 360 - psi_deg + center[1]) ** 2)


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

    # output status report
    s = 'Ran %6d steps [%s]. Found states [%s]' % (
        len(traj) - 1,
        (len(traj) - 1) * engine.snapshot_timestep,
        found_states)
    print s

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
        sample_set=total_sample_set,
        strategies=[
            'extend-complex',
            'extend-minimal'],
        engine=engine)

# -----------------------------------------------------------------------------
# Equilibration
# -----------------------------------------------------------------------------
print """Equilibration"""


# Equilibration
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
