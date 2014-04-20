from distutils.core import setup
from distutils.extension import Extension 
from Cython.Distutils import build_ext
import numpy 

ext_modules = [Extension("Potential_MM_Cython",["Potential_MM_Cython.pyx"], 
    include_dirs = [numpy.get_include()],
    libraries = ["m"] 
    )]

setup(
        name = "Potential_MM_Cython.app",
        cmdclass = {'build_ext':build_ext},
        ext_modules = ext_modules
)
