#!/usr/bin/env python
# create an observation centered at the north galactic pole (sparse!)

import SPyFFI.Observation
o = SPyFFI.Observation.SkyFFI(ra=192.85948, dec=27.12830, label='20150914')
ccds = o.camera.ccds*1
for c in ccds:
    o.camera.ccds = [c]
    o.camera.catalog.addLCs(fmax=0.1, magmax=None)
    o.create()
