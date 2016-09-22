# coding: utf-8# Run from bootstrap paths





import openpathsampling as paths
import numpy as np
import math
### Loading things from storage



old_store = paths.AnalysisStorage("ala_mstis_bootstrap.nc")
print "PathMovers:", len(old_store.pathmovers)
print "Engines:", len(old_store.engines)
print "Samples:", len(old_store.samples)
print "Trajectories:", len(old_store.trajectories)
print "Ensembles:", len(old_store.ensembles)
print "SampleSets:", len(old_store.samplesets)
print "Snapshots:", len(old_store.snapshots)
print "Networks:", len(old_store.networks)
# template = old_store.snapshots[0]
engine = old_store.engines['default']
mstis = old_store.networks[0]
sset = old_store.tag['sampleset']
sset.sanity_check()
## Running RETIS



# logging creates ops_output.log file with details of what the calculation is doing
#import logging.config
#logging.config.fileConfig("logging.conf", disable_existing_loggers=False)
storage = paths.storage.Storage("ala_mstis_production.nc", "w")
storage.snapshots.save(template)
scheme = paths.DefaultScheme(mstis, engine)
mstis_calc = paths.PathSampling(
    storage=storage,
    sample_set=sset,
    move_scheme=scheme
)
mstis_calc.save_frequency = 10
mstis_calc.run(20)
print len(storage.steps)
# commented out during development, so we can "run all" and then do more
storage.close()

