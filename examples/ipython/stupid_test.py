#!/usr/bin/env python
import openpathsampling as paths

cv_x = paths.CV_Function("cv_x", lambda x : x.xyz[0][0])
volume = paths.CVRangeVolume(cv_x, 0.0, 1.0)

print volume.name

import numpy as np
import sys
print [np.__version__, np.__file__, sys.version, sys.executable]

print sys.path
