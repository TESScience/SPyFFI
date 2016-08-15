from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

VERSION = '1.0.4'

class numpy_build_options(build_ext):
    """We can't important numpy directly because it's not ready yet."""
    def finalize_options(self):
        build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())


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
    license='MIT',
    url='https://github.com/TESScience/SPyFFI',
    setup_requires=['numpy'],
    install_requires=['zachopy', 'matplotlib', 'numpy', 'astropy==1.1.2', 'astroquery', 'sh'],
    cmdclass={'build_ext':numpy_build_options},
    ext_modules=[Extension("SPyFFI.cosmical_realistic._cosmical",
                           ["cosmical_realistic/_cosmical.c",
                            "cosmical_realistic/cosmical.c",
                            "cosmical_realistic/twister.c",
                            "cosmical_realistic/seed_tw_ran.c",
			    "cosmical_realistic/fmemopen.c"])],
    # Uncomment this if there's a tagged release that's the same as VERSION, and SPyFFI is publicly released
    download_url = 'https://github.com/TESScience/SPyFFI/tarball/{}'.format(VERSION),
)
