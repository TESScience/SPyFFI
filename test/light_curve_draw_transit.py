#!/usr/bin/env python

# noinspection PyUnresolvedReferences
from SPyFFI.Lightcurve import draw_transit
from numpy.random import RandomState
from pprint import pprint

pprint([[draw_transit(prng=RandomState(seed)).model(i) for seed in range(50)]
        for i in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]])
