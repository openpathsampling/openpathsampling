# coding: utf-8
import openpathsampling as paths


# =============================================================================
# RESTART FROM EXISTING PRODUCTION FILE AND CONTINUE
# =============================================================================
print """RESTART FROM EXISTING PRODUCTION FILE AND CONTINUE"""

print """# Restarting to generate more data"""

storage = paths.Storage("ala_mstis_production.nc", "a")
mstis_calc = paths.PathSampling(
    storage=storage,
    sample_set=storage.steps.last.active,
    move_scheme=storage.schemes[0]
)
mstis_calc.storage.snapshots.only_mention = True
mstis_calc.run(10000)
print 'Final filesize', storage.file_size_str

storage.close()
