#!/usr/bin/env python
import SPyFFI.Observation
ncp = SPyFFI.Observation.SkySubarray(ra=83.0016666667, dec=-0.299095, subarray=170)
ncp.create(cadencestodo={2:20,20:20,120:1,1800:1}, label='withcosmics', correctcosmics=False)
