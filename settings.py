'''global settings required needed for TESS SPyFFI simulations.'''

from imports import *

# okay, so, we need to specify where all the SPyFFI data will be
# this will...
#  ...first try to find an environment variable $SPYFFIDATA
#  ...then default to the current working directory
# if those directories don't contain the required input data, it will complain!

# load the environment variable
prefix = os.getenv('SPYFFIDATA')

# start with the defaults
if prefix is not None:
    print("$SPYFFIDATA is set to {}".format(prefix))
else:
    prefix = os.path.join(os.path.abspath('.'), 'data/')

    print('''
            The environment variable $SPYFFIDATA does not seem to be set;
            defaulting to "data/" inside your current working directory:
            {}

            '''.format(prefix))

    assert(('n' in raw_input('''Please type "n" if that's not okay.''')) == False)
    zachopy.utils.mkdir(prefix)

# create dirs that will store inputs and outputs
dirs = dict(	plots=os.path.join(prefix, 'plots/'),
				inputs=os.path.join(prefix, 'inputs/'),
				outputs=os.path.join(prefix, 'outputs/'),
				intermediates=os.path.join(prefix, 'intermediates/'))

# make sure all those directories exist
for k in dirs.keys():
	zachopy.utils.mkdir(dirs[k])

# shortcuts
plots = dirs['plots']
inputs = dirs['inputs']
outputs = dirs['outputs']
intermediates = dirs['intermediates']
