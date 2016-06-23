#!/usr/bin/env python
# create an observation centered at the north ecliptic pole (midlatitude)
from ..Observation import Observation, default

# start from the default settings
inputs = default

inputs['camera']['dirprefix'] = 'tests/'
inputs['camera']['label'] = 'compresstest'
inputs['catalog']['name'] = 'testpattern'
inputs['camera']['subarray'] = None
inputs['expose']['compress'][2] = True
inputs['expose']['compress'][120] = True
inputs['expose']['skipcosmics'] = True
inputs['observation']['cadencestodo'] = { 2:2, 120:2, 1800:2}
o = Observation(inputs)

o.create()
