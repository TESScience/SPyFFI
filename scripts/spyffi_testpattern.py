#!/usr/bin/env python
import SPyFFI.Observation
<<<<<<< HEAD
ncp = SPyFFI.Observation.TestPattern(subarray=4196, random=True, magnitudes=[10])
ncp.create(todo={2:3,120:3,1800:3}, label='')



=======
ncp = SPyFFI.Observation.TestPattern(subarray=100)
ncp.expose(correctcosmics=False)
>>>>>>> 1143fddf6a0636abc70e91e1f0e73ebcbfb5aaff
