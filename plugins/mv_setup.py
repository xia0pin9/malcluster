from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

sourcefiles=['mvhash.pyx']
ext_modules=[Extension("mvhash", 
	    sourcefiles,
	    include_dirs=[numpy.get_include()]
            )]

setup(
    name = 'mvhash',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
    )
