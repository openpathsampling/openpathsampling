from __future__ import absolute_import
import numpy as np
import random


# Set one rng for all of OPS
if np.version.version < '1.17':
    # Legacy support, remove when numpy 1.16 support is dropped (Py2)
    class RandomState(np.random.RandomState):
        def __init__(self, *args, **kwargs):
            super(RandomState, self).__init__(*args, **kwargs)

        # method overrides to mimic numpy.random.Generator
        def random(self, *args, **kwargs):
            return self.random_sample(*args, **kwargs)

        def integers(self, *args, **kwargs):
            return self.randint(*args, **kwargs)
    rng = RandomState()
else:
    rng = np.random.default_rng()


def default_rng():
    return rng
