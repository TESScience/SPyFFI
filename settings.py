'''global settings required needed for TESS SPyFFI simulations.'''

from imports import *

# define a folder that all data is held
prefix = os.getenv('SPYFFIDATA')

# if that path doesn't exist, complain!
if not os.path.exists(prefix):
	print ("""
			UH-OH! SPyFFI is trying to use '{}' as its base directory.
			However, it seems that directory doesn't exists. Please make
			your SPYFFIDATA environment variable to a real directory,
			or enter the absolute path to such a directory here (although
			you will have to do this every time, so it's probably better
			just to set your environemnt variable).
			""")

	prefix = raw_input('[enter SPyFFI data path]:)
	assert(os.path.exists(prefix))

# create dirs that will store inputs and outputs
dirs = dict(	plots=prefix + 'plots/',
				inputs=prefix + 'inputs/',
				outputs=prefix + 'outputs/',
				intermediates=prefix + 'intermediates/')

# make sure all those directories exist
for k in dirs.keys():
	zachopy.utils.mkdir(dirs[k])
