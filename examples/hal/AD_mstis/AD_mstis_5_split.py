# coding: utf-8
from AD_mstis_S_split import strip_snapshots
import config as cf

# =============================================================================
# RUN FROM BOOTSTRAP PATHS
# =============================================================================
print """RUN FROM BOOTSTRAP PATHS"""

strip_snapshots(cf.storage_production, cf.storage_split)