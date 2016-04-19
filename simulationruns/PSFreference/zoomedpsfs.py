#!/usr/bin/env python
# create an observation centered at the north ecliptic pole (midlatitude)
from SPyFFI.Observation import Observation, default
import numpy as np

# start from the default settings
inputs = default

binby = 2
for binby in [2,4,8]:
    magnification = 101.0/binby
    for focus in inputs['camera']['psfkw']['focus_toinclude']:
        for stellartemp in inputs['camera']['psfkw']['stellartemp_toinclude']:

            inputs['camera']['dirprefix'] = 'PSFreference/'

            inputs['camera']['label'] = 'focus{:.0f}_stellartemp{:.0f}_magnifyby{:.2f}'.format(focus, stellartemp, magnification).replace('.','p')
            inputs['catalog']['name'] = 'testpattern'
            inputs['catalog']['testpatternkw']['spacing'] = 1000.0

            inputs['expose']['skipcosmics'] = True
            inputs['expose']['jitter'] = False
            inputs['camera']['aberrate'] = False
            inputs['camera']['variablefocus'] = True

            inputs['expose']['display'] = False
            inputs['expose']['writenoiseless'] = False
            inputs['expose']['writecosmics'] = False

            inputs['observation']['cadencestodo'] = {1800:1}
            o = Observation(inputs)


            # made a grid in focalxy positions, that land at the centers of CCD pixels
            x, y = np.meshgrid(np.arange(0, 2049, 256), np.arange(0, 2049, 256))

            x = x.flatten().astype(np.float) + (o.camera.ccds[0].center[0] % 1) - 0.5
            y = y.flatten().astype(np.float) + (o.camera.ccds[0].center[1] % 1) - 0.5

            psf = o.camera.psf
            for c in o.camera.ccds:
                o.camera.cartographer.setCCD(c)
                pos = o.camera.cartographer.point(x,y,'ccdxy')
                cat = o.camera.catalog
                o.camera.catalog.ra, o.camera.catalog.dec = pos.celestial.tuple
                o.camera.catalog.tmag = np.ones_like(o.camera.catalog.ra).flatten()*10.0
                o.camera.catalog.pmra = np.zeros_like(o.camera.catalog.ra).flatten()
                o.camera.catalog.pmdec = np.zeros_like(o.camera.catalog.ra).flatten()
                o.camera.catalog.temperature = np.ones_like(o.camera.catalog.ra).flatten()*stellartemp
                o.camera.catalog.addLCs(fractionofstarswithlc=0.0)

                c.expose(advancecounter=False, **inputs['expose'])

                c.header['MAGNIFY'] = ''
                c.header['MAGNOTE'] = ('', 'this image may be a magnitude PSF reference')
                c.header['MAGFACTO'] = (1, 'PSFs are magnified by this factor')
                c.header['MAGFOCUS'] = (focus, 'all at this focus')
                c.header['MAGTEMP'] = (stellartemp, 'all at this stellartemp')
                c.header['IMAGNIFY'] = (False, 'are the PSFs magnified in this image?')


                c.note = 'unmagnified_focus{:.0f}_stellartemp{:.0f}_'.format(focus, stellartemp) + c.fileidentifier
                unmagnifiedfilename = c.directory + c.note + '.fits'
                c.writeToFITS(c.starimage, unmagnifiedfilename, savetype=np.int32)


                #for i in range(len(x)):
                zoomedimage = np.zeros_like(c.image)
                for i in range(len(x)):
                    pos = o.camera.cartographer.point(x[i], y[i], 'ccdxy')
                    zpsf, zx, zy = psf.magnifiedPSF(pos, stellartemp=stellartemp, focus=focus, binby=binby)
                    ok = (zx >= c.xmin)*(zx < c.xmax)*(zy >= c.ymin)*(zy < c.ymax)
                    zoomedimage[zy[ok], zx[ok]] += zpsf[ok]*c.camera.cadence*c.photons(o.camera.catalog.tmag[i])
                    #c.ds9.one(zoomedimage, frame=2)

                c.header['MAGNIFY'] = ''
                c.header['MAGNOTE'] = ('', 'this image may be a magnitude PSF reference')
                c.header['MAGFACTO'] = (np.float(psf.nsubpixelsperpixel)/binby, 'PSFs are magnified by this factor')
                c.header['MAGFOCUS'] = (focus, 'all at this focus')
                c.header['MAGTEMP'] = (stellartemp, 'all at this stellartemp')
                c.header['IMAGNIFY'] = (True, 'are the PSFs magnified in this image?')


                c.note = 'magnified_focus{:.0f}_stellartemp{:.0f}_'.format(focus, stellartemp)  +c.fileidentifier
                magnifiedfilename = c.directory + c.note + '.fits'
                c.writeToFITS(zoomedimage, magnifiedfilename, savetype=np.int32)
