#!/usr/bin/env python
# create an observation centered at the Orion (crowded!)

import SPyFFI.Observation
o = SPyFFI.Observation.SkyFFI(ra=82.0, dec=1.0, label='fordoty')
ccds = o.camera.ccds*1
for c in ccds:
    o.camera.ccds = [c]
    o.camera.catalog.addLCs(fmax=1.0, magmax=None)
    o.create(cadencestodo={2:1})
