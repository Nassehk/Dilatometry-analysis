import sys
if sys.platform[:3]=='win':
    try:
        from setuptools import setup
        from setuptools import Extension
    except ImportError:
        from distutils.core import setup
        from distutils.extension import Extension
else:
    from distutils.core import setup
from Cython.Build import cythonize
from distutils.core import setup
from Cython.Build import cythonize

import numpy
setup(
    ext_modules = cythonize("master_fitter_utils_03.pyx"),
    include_dirs=[numpy.get_include()]
)
