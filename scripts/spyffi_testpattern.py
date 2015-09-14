#!/usr/bin/env python
import SPyFFI.Observation
ncp = SPyFFI.Observation.TestPattern(subarray=4196, random=True, magnitudes=[10])
ncp.create(todo={2:3,120:3,1800:3}, label='')



