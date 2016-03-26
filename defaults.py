# the inputs go into creating the jittering object
jitterkw = dict(

    # if jitterrms is set to a numerical value,
    #   the code will rescale so that sqrt(dx**2 + dy**2) = jitterrms
    # if jitterrms is set to None,
    #   the code will use the input jitter timeseries as is
    jitterrms = None,

    # the code looks for a jitter timeseries (at any cadence faster than 2s)
    #   located in '$SPYFFIDATA/inputs/{rawjitterbasename}'
    rawjitterbasename="AttErrTimeArcsec_80k.dat"

)

# the inputs go into creating the camera
camerakw = dict(

    label='testing',

    # at what cadence (in seconds) does the camera start?
    cadence=1800,

    # what is the commanded central ra and dec of the field?
    ra  =   270.0,
    dec =    66.56070833333332,

    # what is the angle of +y at field center? (measured from N, through E)
    positionangle = 0.0,
    # NOT YET IMPLEMENTED!


    # if subarray = an integer number
    #   create a square subarray, with that many pixels on a side
    # if subarray = None,
    #   don't mess with subarrays, simply creating four separate CCDs
    subarray = None,

    # should this be a camera that magnifies aberration?
    warpspaceandtime=False,

    # by how many steps should the counter advance each exposure
    #  (set to integer > 1 to speed up time)
    counterstep=1,

    # how many fake postage stamps per CCD should be made?
    stamps = None
    # (ultimately, this should be fleshed out into Stamper object, with options)
)

# these inputs create the catalog of target stars
catalogkw = dict(

    # what type of a catalog is this?
    #   options are ['sky', testpattern']
    name = 'sky',

    # the default settings for a real star catalog
    skykw = dict(fast=True),

    # the default settings for a testpattern catalog
    testpatternkw = dict(
                    # how far apart are stars from each other (")
                    spacing=500.0,
                     # list of min, max magnitudes
                    magnitudes=[6,16],
                    # how far to nudge stars (")
                    randomizenudgesby = 21.1,
                    # random prop. mot. (mas/yr)
                    randomizepropermotionsby = 0.0,
                    # randomize the magnitudes?
                    randomizemagnitudes=False,
                    ),


)

# these keywords will be passed to CCD.expose
exposekw = dict(
    # should the exposures be jittered?
    jitter=True,

    # should the exposures write out to file(s)?
    writesimulated=True,

    # should readout smear be included?
    smear=True,

    # what kind of cosmics should be included? (need to be option?)
    cosmicsversion='fancy',

    # should diffusion of cosmics be done?
    diffusion=True,

    # should we skip cosmic injection?
    skipcosmics=False,

    # should we pretend cosmics don't exist?
    correctcosmics=False,

    # should we write an image of the cosmic rays?
    writecosmics=True,

    # should we write an image with no noise?
    writenoiseless=True,

    jitterscale=1.0, # should we rescale the jitter?

    display=True,
)
observationkw = dict(

    # dictionary of cadences to expose (3 each of 2s, 120s, 1800s exposures)
    cadencestodo = {2:3, 120:3, 1800:3},

    # if collate is True,
    #   ccds will expose in order [1,2,3,4,1,2,3,4,....]
    # if collate is False,
    #   ccds will expose in order [1,1,1,...,2,2,2,...,3,3,3,...,4,4,4...]
    collate = True,
    # type
)

inputs = dict(
                observation=observationkw,
                catalog=catalogkw,
                camera=camerakw,
                jitter=jitterkw,
                expose=exposekw
)
