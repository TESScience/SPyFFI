#!/usr/bin/env python
# create an observation centered orion, with two opposite CCDs

import SPyFFI.Observation
o = SPyFFI.Observation.SkyFFI(ra=82.0, dec=1.0, label='SIM009_dvadrift100x', warpspaceandtime=0.01)
o.camera.catalog.addLCs(fmax=1.0, magmax=None, seed=0)

n = 1315
ccds = o.camera.ccds*1
for i in [0,2]:
    c = ccds[i]
    o.camera.ccds = [c]
    o.create(todo={120:n*15}, stamps=500, skipcosmics=True)
    c.aberrator.plotPossibilities()
    o.create(todo={1800:n}, skipcosmics=True)
