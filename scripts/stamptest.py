#!/usr/bin/env python
# create an observation centered orion, with two opposite CCDs

import SPyFFI.Observation
o = SPyFFI.Observation.SkyFFI(ra=82.0, dec=1.0, label='SIM009_dvadrift')
o.camera.catalog.addLCs(fmax=0.0, magmax=None, seed=0)

ccds = o.camera.ccds*1
for i in [0,2]:
    c = ccds[i]
    o.camera.ccds = [c]
    o.create(todo={120:1}, stamps=500, skipcosmics=True)#500)
    c.aberrator.plotPossibilities()
