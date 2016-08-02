from numpy.random import RandomState
from numpy import float32


class DebugRandomState(RandomState):
    """A PseudoRandom Number Generator that makes uniform random numbers that are platform independent.
    Intended to be used for testing."""
    def __init__(self, *args, **kwargs):
        super(DebugRandomState, self).__init__(*args, **kwargs)

    def uniform(self, low=0.0, high=1.0, size=None):
        return float32(super(DebugRandomState, self).uniform(low=low, high=high, size=size))
