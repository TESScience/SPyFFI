# the inputs go into creating the jittering object
jitterkw = dict(

    # if jitterrms is set to a numerical value,
    #   the code will rescale so that sqrt(dx**2 + dy**2) = jitterrms
    # if jitterrms is set to None,
    #   the code will use the input jitter timeseries as is
    jitterrms = None,

    # the code looks for a jitter timeseries (at any cadence faster than 2s)
    #   located in '$SPYFFIDATA/inputs/{rawjitterbasename}'
    rawjitterbasename="cartoon.jitter",

    # by what factor should we rescale the jitter between exposures (without affecting intraexposure jitter)
    amplifyinterexposurejitter=1.0,
)


# the inputs go into creating the PSF library
psfkw = dict(

    # version name for this basic PSF
    version='RRUasbuilt',

    # what prefix was used for
    debprefix = 'woods_prf_feb2016/RAYS_ptSrc_wSi_Oct27model_AP40.6_75C_F3p314adj',

    nsubpixelsperpixel = 101,
    npixels = 21,

    # in the unbinned PSFs
    focus_toinclude = [0,10],
    stellartemp_toinclude = [4350],#[6440, 3850],#[4350],

    npositions_toinclude = 11,#,21,
    noffsets_toinclude = 11

    # in the binned PSFs
    #[6030,4350]# #[7200, 6030, 4350, 3240]
    #7200.,  6440.,  6030.,  5770.,  5250.,  4350.,  3850.,  3240.
    # take every N-
)

# the inputs that go into creating the focus
focuskw = dict(     # what's the range of allowed focus?
                    span=[0.0, 10.0])

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
    warpspaceandtime = False,

    # should this camera change focus over time?
    variablefocus = True,
    # NOT YET IMPLEMENTED!

    # by how many steps should the counter advance each exposure
    #  (set to integer > 1 to speed up time)
    counterstep=1,

    # how many fake postage stamps per CCD should be made, for each cadence?
    # three options:
    #   if None, then create a full-frame image, for that cadence
    #   if an integer, then create a randomize catalog of postage stamps
    #   if a string, then interpret as a filename containing RA and Dec positions (absolute path!)
    stamps = {2:None, 20:None, 120:None, 1800:None},
    
    # (ultimately, this should be fleshed out into Stamper object, with options)

    # include the PSF keywords here, so they can be passed to PSF
    psfkw = psfkw,

    # include the jitter keywords here, so they can be passed to jitter
    jitterkw = jitterkw,

    # include the focus keywords here
    focuskw = focuskw
)

# these inputs create the catalog of target stars
catalogkw = dict(

    # what type of a catalog is this?
    #   options are ['sky', testpattern']
    name = 'sky',

    # should we populate stars with light curves?
    starsarevariable = True,

    # the default settings for a real star catalog
    skykw = dict(fast=False, faintlimit=None),

    # the default settings for a testpattern catalog
    testpatternkw = dict(
                    # how far apart are stars from each other (")
                    spacing=500.0,
                     # list of min, max magnitudes
                    magnitudes=[10,10],
                    # how far to nudge stars (")
                    randomizenudgesby = 21.1,
                    # random prop. mot. (mas/yr)
                    randomizepropermotionsby = 0.0,
                    # randomize the magnitudes?
                    randomizemagnitudes=False,
                    ),

    # the default settings for the light curves for the catalog
    lckw = dict(    # what kinds of variabilty should be allowed?
                    options=['trapezoid', 'sin'],

                    # what magnitude is the faintest star that gets an lc
                    fainteststarwithlc=None,

                    # what fraction of the bright-enough stars get light curves?
                    fractionofstarswithlc=0.5,

                    # (the following keywords get fed into "random()")
                    # what fraction of light curves are extreme?
                    fractionwithextremelc=0.005,

                    # what fraction of light curves get trapezoids (or None, for default)
                    fractionwithtrapezoid=0.3,

                    # what fraction of light curves get sin curves (or None, for default)
                    fractionwithrotation=0.2,

                    # what fraction of light curves get custom light curves (from 0 to 1)
                    fractionwithcustom=0.1,

                    # a seed for the randomizer, for repeatability
                    seed=0)

)

# these keywords will be passed to CCD.expose
exposekw = dict(


    # should the exposures write out to file(s)?
    writesimulated=True,

    # should we write an image of the cosmic rays?
    writecosmics=False,

    # should we write an image with no noise?
    writenoiseless=False,

    # should we compress the images when writing images?
    compress={2:True, 20:True, 120:True, 1800:False},

    # down to what magnitudes should we include? (for fast testing)
    magnitudethreshold=999,

    # should the exposures be jittered?
    jitter=True,

    # should readout smear be included?
    smear=False,

    # should we skip cosmic injection?
    skipcosmics=True,

    # what kind of cosmics should be included? (need to be option?)
    cosmicsversion='fancy',

    # should diffusion of cosmics be done?
    cosmicsdiffusion=True,

    # should we pretend cosmics don't exist?
    correctcosmics=True,

    # should we display images in ds9, as they're created?
    display=False,
)
observationkw = dict(

    # dictionary of cadences to expose (3 each of 2s, 120s, 1800s exposures)
    cadencestodo = {2:3, 20:3, 120:3, 1800:3},

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
