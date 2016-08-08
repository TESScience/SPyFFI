from setuptools import setup, Extension
import numpy.distutils.misc_util

VERSION = '1.0.0'

setup(
    name='SPyFFI',
    version=VERSION,
    package_dir={
        'SPyFFI': '.',
        'SPyFFI.debug': './debug',
        'SPyFFI.documentation': './documentation',
    },
    packages=['SPyFFI', 'SPyFFI.debug', 'SPyFFI.documentation', 'SPyFFI.cosmical_realistic'],
    package_data={
        'SPyFFI': [
            'data/*.txt',  # Light Curve tables
            'relations/*.txt',  # Flux data
        ],
    },
    description="Spiffy Python for Full Frame Images (tools for simulating TESS images, calibrations, and light curves)",
    author='Zach Berta-Thompson',
    author_email='zkbt@mit.edu',
    url='https://github.com/TESScience/SPyFFI',
    install_requires=['zachopy', 'matplotlib', 'numpy', 'astropy==1.1.2', 'astroquery'],
    ext_modules=[Extension("SPyFFI.cosmical_realistic._cosmical",
                           ["cosmical_realistic/_cosmical.c",
                            "cosmical_realistic/cosmical.c",
                            "cosmical_realistic/twister.c",
                            "cosmical_realistic/seed_tw_ran.c",
			    "cosmical_realistic/fmemopen.c"],
                           include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs())],
    # Uncomment this if there's a tagged release that's the same as VERSION, and SPyFFI is publicly released
    download_url = 'https://github.com/TESScience/SPyFFI/tarball/{}'.format(VERSION),
)
