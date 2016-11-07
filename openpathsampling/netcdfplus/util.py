__author__ = 'jan-hendrikprinz'

from time import time as tt

import logging

# logging.basicConfig()
logger = logging.getLogger(__name__)

enable_timing = True


def with_timing_logging(func):
    if enable_timing:
        def _wrapped(*args, **kwargs):
            t1 = tt()
            result = func(*args, **kwargs)
            t2 = tt()
            logger.info('Ran %s in time %f' % (func.func_name, t2 - t1))
            return result

        return _wrapped
    else:
        return func
