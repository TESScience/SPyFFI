#!/usr/bin/env python
import SPyFFI.Observation
ncp = SPyFFI.Observation.SkySubarray(ra=270.0, dec=66.56070833333332)
ncp.camera.catalog.addLCs()
#ncp = SPyFFI.Observation.SkyFFI(ra=82.0, dec=1.0)
ncp.create(todo={1800:100}, label='crowded')
