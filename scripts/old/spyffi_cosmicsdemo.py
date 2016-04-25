#!/usr/bin/env python
import SPyFFI.Observation
ncp = SPyFFI.Observation.SkySubarray(ra=270.0, dec=66.56070833333332, cadence=120)
ncp.expose(correctcosmics=False)
ncp.expose(correctcosmics=True)
