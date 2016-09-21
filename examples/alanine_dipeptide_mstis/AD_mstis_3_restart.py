# coding: utf-8# Restart from existing production file and continue

#### Restarting to generate more data

import openpathsampling as paths
storage = paths.Storage("ala_mstis_production.nc", "a")
mstis_calc = paths.PathSampling(
    storage=storage,
    sample_set=storage.steps.last.active,
    move_scheme=storage.schemes[0]
)
mstis_calc.storage.snapshots.only_mention = True
mstis_calc.run(100)
storage.file_size_str
storage.close()

