#!/usr/bin/env python
from ..Observation import Observation, default

# initialize to default settings
inputs = default

# provide a label, that sets the directory in which outputs will be stored
inputs['camera']['label'] = 'teststamps'


# if the catalog name is set to 'testpattern', create a uniformly space grid
inputs['catalog']['name'] = 'testpattern'


# cadencestodo should be a dictionary of cadences to expose, for example:
# "{2:3, 120:3, 1800:3}" generates (3 each of 2s, 120s, 1800s exposures)
inputs['observation']['cadencestodo'] = {2:3, 120:3, 1800:3}

# (this links closely to cadences to do)
# stamps should be a dictionary, with cadences as keys
# if a cadence's entry is:
#   None        -- a full-frame image will be produced
#   an integer  -- this number of postage stamps will be randomly placed
#   a string    -- this will be interpreted as a filename pointing to a
#                   three-column ascii text file to define where the stamps
#                   should be placed. the columns should be:
#                       [1] RA (in degrees)
#                       [2] Dec (in degrees)
#                       [3] radius (in pixels) of postage stamp
#           the filename path is relative to where you from this script from,
#           or an absolute path. it should include all postage stamps for the
#           entire camera (all four CCDs).
inputs['camera']['stamps'] = {2:'example.stamps', 120:400, 1800:None}


# a dictionary, like those above
# should exposures of a particular cadence be compressed?
inputs['expose']['compress'] = {2:False, 120:False, 1800:False}



'''
------------_--------------_--------------_--------------_--------------_-------
finally, create an observation object, using all these inputs, and make images!
------------_--------------_--------------_--------------_--------------_-------
'''
# generate the observation object
o = Observation(inputs)
# use that object to perform all the exposures
o.create()
