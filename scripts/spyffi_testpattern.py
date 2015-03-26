#!/usr/bin/env python
import SPyFFI.Observation
ncp = SPyFFI.Observation.TestPattern(subarray=100)
ncp.expose(correctcosmics=False)
