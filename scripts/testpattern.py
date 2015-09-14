#!/usr/bin/env python
# create an observation centered at the north ecliptic pole (midlatitude)
import matplotlib
matplotlib.use('agg')
import SPyFFI.Observation
o = SPyFFI.Observation.TestPattern(subarray=4196, random=True, magnitudes=[6,16], label='20150914')
o.camera.catalog.addLCs(fmax=0.1, magmax=None)
o.create()
