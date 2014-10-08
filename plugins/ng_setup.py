from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

sourcefiles=['nghash.pyx']
ext_modules=[Extension("nghash", 
	    	sourcefiles
# 	    	extra_compile_args=["-O3", "-march=native", "-ffast-math","-funroll-loops","-D NPY_NO_DEPRECATED_API"],
# 	    	define_macros=[("NPY_NO_DEPRECATED_API", None)],
#  	    	include_dirs = [numpy.get_include()]
            )]

setup(
    name = 'nghash',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [numpy.get_include()],
    ext_modules = ext_modules
    )
