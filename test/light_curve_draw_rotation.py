#!/usr/bin/env python

from __future__ import print_function

# noinspection PyUnresolvedReferences
from SPyFFI.Lightcurve import draw_rotation
from numpy.random import RandomState
import json

print(json.dumps(
    [[draw_rotation(prng=RandomState(seed)).model(i) for seed in range(10)]
     for i in [0.0, 0.1, 0.2, 0.3, 0.4]],
    indent=4))
