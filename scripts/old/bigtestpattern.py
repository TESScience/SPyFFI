#!/usr/bin/env python
# create an observation centered at the north ecliptic pole (midlatitude)

import SPyFFI.Observation
o = SPyFFI.Observation.TestPattern(subarray=4196, random=True, magnitudes=[6,16], label='20150914')
o.camera.catalog.addLCs(fmax=1.0, magmax=None)
o.create()
