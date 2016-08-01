from setuptools import setup, Extension
import numpy.distutils.misc_util

setup(install_requires=['numpy'],
      ext_modules=[Extension("_cosmical", ["_cosmical.c", "cosmical.c", "twister.c", "seed_tw_ran.c"],
                             include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs())], )
