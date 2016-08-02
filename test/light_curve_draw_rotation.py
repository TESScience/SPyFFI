#!/usr/bin/env python

# noinspection PyUnresolvedReferences
from SPyFFI.Lightcurve import draw_rotation
# noinspection PyUnresolvedReferences
from SPyFFI.debug.rng import DebugRandomState
from pprint import pprint


pprint([[draw_rotation(prng=DebugRandomState(seed)).model(i) for seed in range(10)]
        for i in [0.0, 0.1, 0.2, 0.3, 0.4]])
