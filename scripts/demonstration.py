#!/usr/bin/env python
from ..Observation import Observation, default

# initialize to default settings
inputs = default


'''
------------_--------------_--------------_--------------_--------------_-------
"camera" inputs change many of the basics of the observations.
------------_--------------_--------------_--------------_--------------_-------
'''


# provide a label, that sets the directory in which outputs will be stored
inputs['camera']['label'] = 'demonstration'

# what is the commanded central ra and dec of the field?
inputs['camera']['ra'] = 82.0
inputs['camera']['dec'] = 1.0

# if subarray = an integer number
#   create a square subarray, with that many pixels on a side
# if subarray = None,
#   simply creating four separate CCDs, with their default sizes
inputs['camera']['subarray'] = None



'''
------------_--------------_--------------_--------------_--------------_-------
"catalog" inputs affect what stars will be used to populate the images.
------------_--------------_--------------_--------------_--------------_-------
'''
# if the catalog name is set to 'sky', draw stars from the real sky (UCAC4)
#inputs['catalog']['name'] = 'sky'
#inputs['catalog']['skykw']['faintlimit'] = 10.0

# if the catalog name is set to 'testpattern', create a uniformly space grid
inputs['catalog']['name'] = 'testpattern'


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
inputs['jitter']['jitterrms'] = None

# this will amplify the jitter between exposures (without reblurring the PSFs)
inputs['jitter']['amplifyinterexposurejitter'] = 1.0

# differential velocity aberration will be included. if you want to make
# the smooth positional trends associated with DVA larger, then you can
# set warpspaceandtime=0.1 to slow down the speed of light by a factor of 10X,
# thereby making the aberration 10X larger
# if warpspaceandtime is set to False, the aberration will be as expected
inputs['camera']['warpspaceandtime'] = False

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
inputs['expose']['jitter'] = True

# should readout smear be included?
inputs['expose']['smear'] = False

# should we skip cosmic injection?
inputs['expose']['skipcosmics'] = True

# should we pretend cosmics don't exist?
# (if skipcosmics is false and correctcosmics is true,
#  a cosmic ray image will be made but not added to the final image)
inputs['expose']['correctcosmics'] = True

# should we display images in ds9, as they're created?
inputs['expose']['display'] = False

# should we compress certain image sizes?
inputs['expose']['compress']={2:True, 120:False, 1800:False, 20:True}

'''
------------_--------------_--------------_--------------_--------------_-------
"observation" keywords set the overall group of exposures to be made.
------------_--------------_--------------_--------------_--------------_-------
'''

# cadencestodo should be a dictionary of cadences to expose, for example:
# "{2:3, 120:3, 1800:3}" generates (3 each of 2s, 120s, 1800s exposures)
inputs['observation']['cadencestodo'] = {1800:3, 2:3, 120:3}

# (this links closely to cadences to do)
# stamps should be a dictionary, with cadences as keys
# if a cadence's entry is:
#   None        -- a full-frame image will be produced
#   an integer  -- this number of postage stamps will be randomly places
#   a string    -- this will be interpreted as a filename pointing to a
#                   three-column ascii text file to define where the stamps
#                   should be placed. the columns should be:
#                       [1] RA (in degrees)
#                       [2] Dec (in degrees)
#                       [3] radius (in pixels) of postage stamp
inputs['camera']['stamps'] = {2:None, 120:20, 1800:None}

# a dictionary, like those above
# should exposures of a particular cadence be compressed?
inputs['expose']['compress'] = {2:True, 120:True, 1800:False}




'''
------------_--------------_--------------_--------------_--------------_-------
finally, create an observation object, using all these inputs, and make images!
------------_--------------_--------------_--------------_--------------_-------
'''
# generate the observation object
o = Observation(inputs)
# use that object to perform all the exposures
o.create()
