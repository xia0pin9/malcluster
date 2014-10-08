from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

sourcefiles=['mvhash.pyx']
ext_modules=[Extension("mvhash", 
	    sourcefiles
            )]

setup(
    name = 'mvhash',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
    )
