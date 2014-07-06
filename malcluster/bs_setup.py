from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

sourcefiles=['bshash.pyx']
ext_modules=[Extension("bshash", 
	    sourcefiles
            )]

setup(
    name = 'bshash',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
    )
