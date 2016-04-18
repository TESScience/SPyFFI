#!/usr/bin/env python
# create an observation centered at the north ecliptic pole (midlatitude)

import SPyFFI.Observation
o = SPyFFI.Observation.SkySubarray(ra=50, dec=20, label='aberrationtest', subarray=400, counterstep=350, warpspaceandtime=True)
o = SPyFFI.Observation.SkyFFI(ra=50, dec=20, label='aberrationtest')
#o = SPyFFI.Observation.TestPattern(subarray=200, spacing=500.0, label='aberrationtest', ra=270.0, dec=66.0, warpspaceandtime=False, counterstep=3500)
ccds = o.camera.ccds*1
o.camera.catalog.addLCs(fmax=0.0, magmax=None, seed=0)
for i, c in enumerate(ccds):
    o.camera.ccds = [c]
    o.create(cadencestodo={1800:1}, stamps=100, skipcosmics=True)
    c.aberrator.plotPossibilities()

'''
ra,dec,mag,teff = o.camera.catalog.snapshot(o.ccd.bjd)
stars = o.camera.cartographer.point(ra,dec,type='celestial')

dx,dy,delon = [],[], []
import numpy as np
bjds = np.linspace(0,365,1000)+o.ccd.bjd0
for bjd in bjds:
    nudges = o.ccd.aberrations(stars,bjd)
    dx.append(nudges[0])
    dy.append(nudges[1])
    delon.append(o.ccd.delon)

import matplotlib.pyplot as plt
plt.figure()
gs = plt.matplotlib.gridspec.GridSpec(2,2, left=0.15)
bjds -= min(bjds)
# top row is uncorrected
ax = plt.subplot(gs[0,0])
plt.axvline(27.4, color='gray', alpha=1)
ax.plot(bjds,dx, alpha=0.3)
plt.ylabel('velocity\naberration (pixels)')
plt.title('x')
ax = plt.subplot(gs[0,1],sharey=ax,sharex=ax)
ax.plot(bjds,dy, alpha=0.3)
plt.axvline(27.4, color='gray', alpha=1)
plt.xlim(-1+ min(bjds), max(bjds)+1)
plt.title('y')

ax = plt.subplot(gs[1,0])
ax.plot(bjds,dx-np.mean(dx,1).reshape(len(bjds),1), alpha=0.3)
plt.axvline(27.4, color='gray', alpha=1)
plt.xlabel('Time (days)')
plt.ylabel('differential velocity\naberration (pixels)')

ax = plt.subplot(gs[1,1],sharey=ax,sharex=ax)
ax.plot(bjds,dy-np.mean(dy,1).reshape(len(bjds),1), alpha=0.3)
plt.axvline(27.4, color='gray', alpha=1)
plt.xlim(-1 + min(bjds), max(bjds)+1)
plt.xlabel('Time (days)')
plt.savefig(o.ccd.directory+'aberrationoveroneyear.pdf')
'''
