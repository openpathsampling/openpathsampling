import numpy as np


# Set one rng for all of OPS
if np.version.version < '1.17':
    # Legacy support, remove when numpy 1.16 support is dropped (Py2)
    class RandomState(np.random.RandomState):
        # method overrides to mimic numpy.random.Generator
        def random(self, *args, **kwargs):
            return self.random_sample(*args, **kwargs)

        def integers(self, *args, **kwargs):
            return self.randint(*args, **kwargs)
    DEFAULT_RNG = RandomState()
else:
    DEFAULT_RNG = np.random.default_rng()


def default_rng():
    return DEFAULT_RNG
