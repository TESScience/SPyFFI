#!/usr/bin/env python
# create an observation centered at the north ecliptic pole (midlatitude)

import SPyFFI.Observation
o = SPyFFI.Observation.SkyFFI(ra=270.00000, dec=66.56071, label='20150927')
ccds = o.camera.ccds*1
for c in ccds:
    o.camera.ccds = [c]
    o.camera.catalog.addLCs(fmax=0.1, magmax=None)
    o.create()
