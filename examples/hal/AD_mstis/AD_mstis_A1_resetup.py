# coding: utf-8
import openpathsampling as paths

from openpathsampling.tools import refresh_output

import config as cf


# =============================================================================
# RUN FROM BOOTSTRAP PATHS
# =============================================================================
print """RECREATE RUN FROM THE INITIAL ONE"""

# -----------------------------------------------------------------------------
# Loading things from storage
# -----------------------------------------------------------------------------
print """Loading things from storage"""

old_store = paths.Storage(cf.storage_setup)

template = old_store.tag['template']
engine_high = old_store.engines['high']
engine_low = old_store.engines['default']
mstis = old_store.networks[0]

initial_snapshot = old_store.snapshots.first

# sset = old_store.tag['sampleset']
# sset.sanity_check()

# initialize engine
# if we do not select a platform the fastest possible will be chosen
# but we explicitly request to use the one in the config file
platform = cf.platform
engine_high.initialize(platform)

print 'Engine uses platform `%s`' % engine_high.platform

# Running RETIS


# -----------------------------------------------------------------------------
# Create the storage
# -----------------------------------------------------------------------------
print """Create the storage"""

storage = paths.Storage(cf.storage_resetup, 'w')
storage.save(template)
storage.save(engine_low)
storage.save(engine_high)

# -----------------------------------------------------------------------------
# Copy tags
# -----------------------------------------------------------------------------
print """Copy tags"""

for key, value in old_store.tag.iteritems():
    if key != 'sampleset':
        storage.tag[key] = value

# -----------------------------------------------------------------------------
# Set up the `MoveScheme`
# -----------------------------------------------------------------------------
print """Set up the `MoveScheme`"""

scheme = old_store.schemes[0]

# -----------------------------------------------------------------------------
# Initial trajectories
# -----------------------------------------------------------------------------
print """Initial trajectories"""

engine_high.initialize(platform)

print 'High-Engine uses'
print 'platform `%s`' % engine_high.platform
print 'temperature `%6.2f K' % (
    float(engine_high.integrator.getGlobalVariableByName('kT')) / 0.0083144621)
print

states = old_store.tag['states']
vol_state = {state: old_store.volumes[state] for state in states}

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
