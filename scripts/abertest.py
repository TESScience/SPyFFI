#!/usr/bin/env python
# create an observation centered at the north ecliptic pole (midlatitude)

import SPyFFI.Observation
o = SPyFFI.Observation.SkySubarray(ra=0, dec=0, label='aberrationtest', subarray=200)
#o = SPyFFI.Observation.TestPattern(subarray=4000, spacing=2000.0, label='aberrationtest')
ccds = o.camera.ccds*1
for c in ccds:
    o.camera.ccds = [c]
    o.camera.catalog.addLCs(fmax=0.1, magmax=None)
    o.create(todo={1800:1})



ra,dec,mag,teff = o.camera.catalog.snapshot(o.ccd.bjd)
stars = o.camera.cartographer.point(ra,dec,type='celestial')

dx,dy,delon = [],[], []
import numpy as np
bjds = np.linspace(0,365*2,1000)+o.ccd.bjd0
for bjd in bjds:
    nudges = o.ccd.aberrations(stars,bjd)
    dx.append(nudges[0])
    dy.append(nudges[1])
    delon.append(o.ccd.delon)

import matplotlib.pyplot as plt
trend = np.mean(dx,1)
plt.figure()
plt.plot(bjds,dx, alpha=0.3)
plt.axvline(o.ccd.bjd0+27.4)

plt.figure()
plt.plot(bjds,dx-trend.reshape(len(bjds),1), alpha=0.3)
plt.axvline(o.ccd.bjd0+27.4)
