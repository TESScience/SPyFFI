#!/usr/bin/env python
# create an observation centered at the north ecliptic pole (midlatitude)
from ..Observation import Observation, default

# start from the default settings
inputs = default

inputs['camera']['label'] = 'soren'
inputs['catalog']['name'] = 'testpattern'
inputs['camera']['subarray'] = None
inputs['camera']['variablefocus'] = False
inputs['expose']['jitterscale'] = 1.0
inputs['expose']['skipcosmics'] = True
inputs['catalog']['testpatternkw']['randomizemagnitudes'] = True
inputs['observation']['cadencestodo'] = {1800:9}#, 120:16, 1800:16}
o = Observation(inputs)

o.create()
