from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

sourcefiles=['sdhash.pyx']
ext_modules=[Extension("sdhash", 
	    sourcefiles
            )]

setup(
    name = 'sdhash',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
    )
