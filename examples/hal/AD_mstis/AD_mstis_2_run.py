# coding: utf-8
import openpathsampling as paths

import config as cf


# =============================================================================
# RUN FROM BOOTSTRAP PATHS
# =============================================================================
print """RUN FROM BOOTSTRAP PATHS"""

# -----------------------------------------------------------------------------
# Loading things from storage
# -----------------------------------------------------------------------------
print """Loading things from storage"""

old_store = paths.AnalysisStorage(cf.storage_setup)

# template = old_store.snapshots[0]
engine = old_store.engines['default']
mstis = old_store.networks[0]
sset = old_store.tag['sampleset']
sset.sanity_check()

# initialize engine
# if we do not select a platform the fastest possible will be chosen
# but we explicitly request to use the one in the config file
platform = cf.platform
engine.initialize(platform)

print 'Engine uses platform `%s`' % engine.platform

# Running RETIS
storage = paths.storage.Storage("ala_mstis_production.nc", "w")
storage.save(old_store.snapshots[0])

scheme = paths.DefaultScheme(mstis, engine)
mstis_calc = paths.PathSampling(
    storage=storage,
    sample_set=sset,
    move_scheme=scheme
)

mstis_calc.save_frequency = 50
mstis_calc.run(10000)
print 'steps:', len(storage.steps)
print 'snapshots:', len(storage.snapshots)
print 'trajectories:', len(storage.trajectories)
print 'samples:', len(storage.samples)
print 'filesize:', storage.file_size_str

storage.close()
