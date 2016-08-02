"""Global settings required needed for TESS SPyFFI simulations."""
import os
import logging

# okay, so, we need to specify where all the SPyFFI data will be
# this will...
#  ...first try to find an environment variable $SPYFFIDATA
#  ...then default to "~/.tess/spyffi"
# if those directories don't contain the required input data, it will complain!

# load the environment variable
prefix = os.getenv('SPYFFIDATA', os.path.expanduser("~/.tess/spyffi"))
if not os.path.exists(prefix):
    os.makedirs(prefix)


logging.basicConfig(
    format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s",
    datefmt="%H:%M:%S",
    level=getattr(logging, os.getenv('LOG', 'WARNING').upper()))

log_file_handler = logging.FileHandler(os.getenv('LOG_FILE', os.path.join(prefix, 'SPyFFI.log')))

# create dirs that will store inputs and outputs
dirs = {'plots': os.path.join(prefix, 'plots/'), 'inputs': os.path.join(prefix, 'inputs/'),
        'outputs': os.path.join(prefix, 'outputs/'), 'intermediates': os.path.join(prefix, 'intermediates/')}

# make sure all those directories exist
for d in dirs.values():
    if not os.path.exists(d):
        os.makedirs(d)

# shortcuts
plots = dirs['plots']
inputs = dirs['inputs']
outputs = dirs['outputs']
intermediates = dirs['intermediates']
