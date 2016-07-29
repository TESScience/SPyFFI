#!/usr/bin/env python
# create an observation centered at the north ecliptic pole (midlatitude)
from SPyFFI.Lightcurve import draw_transit
from pprint import pprint

# start from the default settings
pprint([[draw_transit(seed=seed).model(i) for seed in range(50)]
        for i in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]])
