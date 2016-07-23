#!/usr/bin/env python
# create an observation centered at the north ecliptic pole (midlatitude)
from SPyFFI.Observation import Observation, default

# start from the default settings
inputs = default

inputs['camera']['label'] = 'smallone'
inputs['camera']['subarray'] = 400
inputs['observation']['testpattern'] = True
inputs['catalog']['name'] = 'sky'
inputs['expose']['skipcosmics'] = True
inputs['expose']['writenoiseless'] = True
inputs['observation']['cadencestodo'] = {1800:1}
o = Observation(inputs)

o.create()
