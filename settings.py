'''Keeps track of settings required needed for TESS simulations, particularly dirs.'''
from imports import *

# define a folder that all data is held
prefix = "/Users/zkbt/Cosmos/Data/TESS/FFIs/"
# if that path doesn't exist, complain!
if not os.path.exists(prefix):
	print "Please enter a parent directory where you want to store all data related to TESS FFI's:"
	a = raw_input('   ')
	if os.path.exists(prefix):
		prefix = a
	else:
		assert(False)

# create dirs that will store inputs and outputs
dirs = dict(plots=prefix + 'plots/', inputs=prefix + 'inputs/', outputs=prefix + 'outputs/', intermediates=prefix + 'intermediates/')
for k in dirs.keys():
	zachopy.utils.mkdir(dirs[k])
