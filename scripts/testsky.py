#!/usr/bin/env python
# create an observation centered at the north ecliptic pole (midlatitude)
from ..Observation import Observation, default

# start from the default settings
inputs = default

inputs['camera']['label'] = 'newtest'
inputs['catalog']['name'] = 'sky'
inputs['expose']['jitterscale'] = 100.0
inputs['expose']['skipcosmics'] = True
inputs['observation']['cadencestodo'] = {1800: 9}  # , 120:16, 1800:16}
o = Observation(inputs)

# o.create()
