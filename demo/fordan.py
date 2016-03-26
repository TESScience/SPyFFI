import SPyFFI.Observation

inputs = dict(subarray=200, label='dfm_extreme',
                ra=82.0, dec=1.0,
                jitter=True,
                warpspaceandtime=0.01, jitterscale=10.0,
                cadence=1800, correctcosmics=True,
                write=True, writenoiseless=False)

o = SPyFFI.Observation.SkySubarray(**inputs)
ccds = o.camera.ccds*1
for c in ccds:
    o.camera.ccds = [c]
    o.camera.catalog.addLCs(fmax=0.1, magmax=None)
    o.create(cadencestodo={1800:1320})
