from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

sourcefiles=['imphash.pyx']
ext_modules=[Extension("imphash", sourcefiles)]

setup(
    name = 'imphash',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [numpy.get_include()],
    ext_modules = ext_modules
    )
