#!/usr/bin/env python
import SPyFFI.Observation
ncp = SPyFFI.Observation.SkyFFI(ra=82.0, dec=1.0)
ncp.create(todo={2:3,120:3,1800:100}, label='crowded')
