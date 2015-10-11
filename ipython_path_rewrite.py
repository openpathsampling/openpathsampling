#!/usr/bin/env python
#minimal script to show how IPython rewrites paths in an environment

import sys
from IPython import start_ipython

print sys.path
start_ipython(['-c', '"import sys; print sys.path"'])
