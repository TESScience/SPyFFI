#!/usr/bin/env python
import SPyFFI.Observation
ncp = SPyFFI.Observation.SkyFF(ra=82.0, dec=1.0)
ncp.create(todo={2:3,120:3,1800:3}, label='crowded')
