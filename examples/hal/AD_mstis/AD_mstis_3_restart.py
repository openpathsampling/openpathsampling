# coding: utf-8
import openpathsampling as paths

import config as cf

# =============================================================================
# RESTART FROM EXISTING PRODUCTION FILE AND CONTINUE
# =============================================================================
print """RESTART FROM EXISTING PRODUCTION FILE AND CONTINUE"""

print """# Restarting to generate more data"""

# -----------------------------------------------------------------------------
# Loading things from storage
# -----------------------------------------------------------------------------
print """Open existing storage"""

storage = paths.Storage(cf.storage_production, "a")

print """Initialize engine to use given platform"""
engine = storage.engines['default']
engine.initialize(cf.platform)
print 'Engine uses platform `%s`' % engine.platform

last_sampleset = storage.steps.last.active
scheme = storage.schemes[0]

mstis_calc = paths.PathSampling(
    storage=storage,
    sample_set=last_sampleset,
    move_scheme=scheme
)
mstis_calc.storage.snapshots.only_mention = True
mstis_calc.run(10000)

print 'steps:', len(storage.steps)
print 'snapshots:', len(storage.snapshots)
print 'trajectories:', len(storage.trajectories)
print 'samples:', len(storage.samples)
print 'filesize:', storage.file_size_str

storage.close()
