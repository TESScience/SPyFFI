"""Global settings required needed for TESS SPyFFI simulations."""
import os
import logging
# noinspection PyUnresolvedReferences
from sh import wget, tar, rm, shasum

# okay, so, we need to specify where all the SPyFFI data will be
# this will...
#  ...first try to find an environment variable $SPYFFIDATA
#  ...then default to "~/.tess/spyffi"
# if those directories don't contain the required input data, it will complain!

# load the environment variable
prefix = os.path.expanduser(os.getenv('SPYFFIDATA', "~/.tess/spyffi"))
if not os.path.exists(prefix):
    os.makedirs(prefix)

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s",
    datefmt="%H:%M:%S",
    level=getattr(logging, os.getenv('LOG', 'WARNING').upper()))

log_file_handler = logging.FileHandler(os.getenv('LOG_FILE', os.path.join(prefix, 'SPyFFI.log')))
logger = logging.getLogger(__name__)
logger.addHandler(log_file_handler)

# create dirs that will store inputs and outputs
dirs = {'plots': os.path.join(prefix, 'plots/'),
        'inputs': os.path.join(prefix, 'inputs/'),
        'outputs': os.path.join(prefix, 'outputs/'),
        'intermediates': os.path.join(prefix, 'intermediates/')}

# TODO: https://github.com/TESScience/SPyFFI/issues/18
checksums = \
    """
994ddb7fb44efd7963730d29051378e0c79861ccd5a93de0bb46cfe34d2b8b8f  ./inputs/AttErrTimeArcsec_80k.dat
454a4219d5d97b53416754c54e37a7e98c071dab09ab1151084ab1d1f0891d95  ./inputs/cartoon.jitter
dc518dd3c58f01f1a723b21ed3ce56bb48b56ae8a99614024a8db5535e3809a4  ./inputs/covey_pickles.txt
494d685fce7d322f9c310f9995207d8977e84a7cfcf0dcab444ba62455980819  ./inputs/davenport_table1.txt
a659f51a89b84b2cb194e83984936fc63f83e7d52f1ff9f3317321a5021e6478  ./inputs/pickles_table2.txt
6cf1fbbff6b604471a947039bc0d7f56b06c69cd234c91cc72bd989225d359bb  ./intermediates/psfs/RRUasbuilt/focus0and10_stellartemp4350/originaldeblibrary.npy
dbdac6e3268267eb286fb49fbccec3d4ba0e34c6a611ea18de4d54961f910cdd  ./intermediates/psfs/RRUasbuilt/focus0and10_stellartemp4350/pixelizedlibrary_cartoon.jitter.cadence120s.unscaled_perfectpixels_11positions_11offsets.npy
dbdac6e3268267eb286fb49fbccec3d4ba0e34c6a611ea18de4d54961f910cdd  ./intermediates/psfs/RRUasbuilt/focus0and10_stellartemp4350/pixelizedlibrary_cartoon.jitter.cadence1800s.unscaled_perfectpixels_11positions_11offsets.npy
dbdac6e3268267eb286fb49fbccec3d4ba0e34c6a611ea18de4d54961f910cdd  ./intermediates/psfs/RRUasbuilt/focus0and10_stellartemp4350/pixelizedlibrary_cartoon.jitter.cadence2s.unscaled_perfectpixels_11positions_11offsets.npy
"""
data_url = "https://www.dropbox.com/s/0e4c2uk34phv4qx/SPyFFI_coreinputs.tar.gz"


def initialize():
    # noinspection PyUnresolvedReferences
    from sh import wget, tar, rm, shasum
    if not os.path.exists(prefix):
        os.makedirs(prefix)

    if (not os.path.exists(dirs['inputs'])) or (not os.path.exists(dirs['intermediates'])):
        try:
            if not os.path.exists(prefix):
                logger.info("Creating {DIR}".format(DIR=prefix))
                os.makedirs(prefix)
            logger.info("Downloading data from {URL} to {DIR}".format(URL=data_url, DIR=prefix))
            tar(wget(data_url, "-qO-", _piped=True), "xz", _cwd=prefix)
            logger.info("Checking checksums of downloaded files")
            for line in shasum("-c", _cwd=prefix, _in=checksums, _iter=True):
                logger.info(line)
        except Exception as e:
            logger.info("Error: {}".format(e.message))
            logger.info("Deleting {DIR}".format(DIR=dirs['inputs']))
            rm(dirs['inputs'], '-rf')
            logger.info("Deleting {DIR}".format(DIR=dirs['intermediates']))
            rm(dirs['intermediates'], '-rf')
            raise

    # make sure all those directories exist
    for d in (dirs['outputs'], dirs['plots']):
        if not os.path.exists(d):
            logger.info("Creating {DIR}".format(DIR=d))
            os.makedirs(d)


initialize()

# shortcuts
plots = dirs['plots']
inputs = dirs['inputs']
outputs = dirs['outputs']
intermediates = dirs['intermediates']
