#!/usr/bin/env python

from __future__ import print_function

# noinspection PyUnresolvedReferences
from SPyFFI.Lightcurve import draw_transit
from numpy.random import RandomState
import json

print(json.dumps([[draw_transit(prng=RandomState(seed)).model(i).tolist() for seed in range(50)]
                  for i in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]],
                 indent=4))
