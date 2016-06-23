from distutils.core import setup, Extension
import numpy.distutils.misc_util

VERSION = '0.1.0'

setup(
    name='SPyFFI',
    version=VERSION,
    package_dir={'spyffi': '.'},
    packages=['spyffi'],
    description="Spiffy Python for Full Frame Images (tools for simulating TESS images, calibrations, and light curves)",
    author='Zach Berta-Thompson',
    author_email='zkbt@mit.edu',
    url='https://github.com/zkbt/SPyFFI',
    # TODO: Version bump this when zachopy has a new release
    install_requires=['zachopy==0.1.0', 'matplotlib', 'numpy', 'astropy', 'astroquery'],
    dependency_links=['git+https://github.com/zkbt/zachopy.git#egg=zachopy-0.1.0'],
    ext_modules=[Extension("cosmical_realistic._cosmical", ["cosmical_realistic/_cosmical.c",
                                                            "cosmical_realistic/cosmical.c",
                                                            "cosmical_realistic/twister.c",
                                                            "cosmical_realistic/seed_tw_ran.c"],
                           include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs())],
    # Uncomment this if there's a tagged release that's the same as VERSION, and SPyFFI is publicly released
    # download_url = 'https://github.com/zkbt/SPyFFI/tarball/{}'.format(VERSION),
)
