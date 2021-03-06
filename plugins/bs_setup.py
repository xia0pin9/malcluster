from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

sourcefiles=['bshash.pyx']
ext_modules=[Extension("bshash", 
	    sourcefiles
            )]

setup(
    name = 'bshash',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [numpy.get_include()],
    ext_modules = ext_modules
    )
