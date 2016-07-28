#!/usr/bin/env python
# create an observation centered at the north ecliptic pole (midlatitude)
from SPyFFI.Lightcurve import draw_rotation
from pprint import pprint

# start from the default settings
pprint([[draw_rotation(seed=seed).model(i) for seed in range(10)]
        for i in [0.0, 0.1, 0.2, 0.3, 0.4]])
