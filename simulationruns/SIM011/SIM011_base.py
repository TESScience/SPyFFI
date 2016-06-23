#!/usr/bin/env python
from ...Observation import default

# initialize to default settings
inputs = default


'''
------------_--------------_--------------_--------------_--------------_-------
"camera" inputs change many of the basics of the observations.
------------_--------------_--------------_--------------_--------------_-------
'''


inputs['camera']['dirprefix'] = 'SIM011and12/'

# provide a label, that sets the directory in which outputs will be stored
inputs['camera']['label'] = 'SIM011and12'

# if subarray = an integer number
#   create a square subarray, with that many pixels on a side
# if subarray = None,
#   simply creating four separate CCDs, with their default sizes
inputs['camera']['subarray'] = None
inputs['camera']['counterstep'] = 1

inputs['camera']['aberrate'] = False
inputs['camera']['warpspaceandtime'] = 0.1

inputs['camera']['variablefocus'] = False


inputs['camera']['psfkw']['npositions_toinclude'] = 11
inputs['camera']['psfkw']['noffsets_toinclude'] = 11

'''
------------_--------------_--------------_--------------_--------------_-------
"catalog" inputs affect what stars will be used to populate the images.
------------_--------------_--------------_--------------_--------------_-------
'''
# if the catalog name is set to 'sky', draw stars from the real sky (UCAC4)
inputs['catalog']['name'] = 'sky'
inputs['catalog']['skykw']['faintlimit'] = None

# if the catalog name is set to 'testpattern', create a uniformly space grid
#inputs['catalog']['name'] = 'testpattern'


'''
------------_--------------_--------------_--------------_--------------_-------
"jitter" inputs change the jitter, with an effect on both exposure-to-exposure
nudges and intra-exposure blurring. They may require (expsensive) recomputation
of the PSF library.
------------_--------------_--------------_--------------_--------------_-------
'''
# the code looks for a jitter timeseries (sampled at any cadence faster than 2s)
#   located in '$SPYFFIDATA/inputs/{rawjitterbasename}'
# inputs['jitter']['rawjitterbasename'] = "AttErrTimeArcsec_80k.dat"

# if jitterrms is set to a numerical value,
#   the code will rescale so that sqrt(dx**2 + dy**2) = jitterrms (in arcsec)
# if jitterrms is set to None,
#   the code will use the input jitter timeseries as is
inputs['jitter']['jitterrms'] = None#2.0/3.0

inputs['jitter']['amplifyinterexposurejitter'] = 10.0

'''
------------_--------------_--------------_--------------_--------------_-------
"expose" keywords determine how individual exposures are generated.
------------_--------------_--------------_--------------_--------------_-------
'''

# should the exposures write out to file(s)?
inputs['expose']['writesimulated'] = True

# should we write an image of the cosmic rays?
inputs['expose']['writecosmics'] = False

# should we write an image with no noise?
inputs['expose']['writenoiseless'] = True

# down to what magnitudes should we include? (for fast testing)
inputs['expose']['magnitudethreshold'] = 999

# should the exposures be jittered?
inputs['expose']['jitter'] = False

# should readout smear be included?
inputs['expose']['smear'] = False

# should we skip cosmic injection?
inputs['expose']['skipcosmics'] = True

# should we pretend cosmics don't exist?
inputs['expose']['correctcosmics'] = True

# should we display images in ds9, as they're created?
inputs['expose']['display'] = False

'''
------------_--------------_--------------_--------------_--------------_-------
"observation" keywords set the overall group of exposures to be made.
------------_--------------_--------------_--------------_--------------_-------
'''

# cadencestodo should be a dictionary of cadences to expose, for example:
# "{2:3, 120:3, 1800:3}" generates (3 each of 2s, 120s, 1800s exposures)
inputs['observation']['cadencestodo'] = {1800:int(24*2*27.6)}#int(24*2*13.7)
inputs['observation']['collate'] = True
