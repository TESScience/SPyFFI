#!/usr/bin/env python
# create an observation centered at the north ecliptic pole (midlatitude)
from ..Observation import Observation, default

# start from the default settings
inputs = default

inputs['camera']['label'] = 'newtest'
inputs['catalog']['name'] = 'testpattern'
inputs['camera']['subarray'] = 100
inputs['expose']['skipcosmics'] = True
inputs['catalog']['testpatternkw']['randomizemagnitudes'] = True
inputs['observation']['cadencestodo'] = {2:16, 120:16, 1800:16}
o = Observation(inputs)


psf = o.camera.psf
'''
for i in np.arange(-2000,2000,100):
    pos = psf.camera.cartographer.point(i, 500, 'focalxy')
    psf.display.one(psf.highResolutionPSF(pos, chatty=True).T, clobber=False)
'''
pos = psf.camera.cartographer.point(300, 1200, 'focalxy')
psf.populateBinned(plot=True)
b = psf.binHighResolutionPSF(pos, chatty=True, plot=True)
