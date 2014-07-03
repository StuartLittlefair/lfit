from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy 
import os

library_dirs = []
include_dirs = []

# need to direct to where includes and  libraries are
if os.environ.has_key('TRM_SOFTWARE'):
    library_dirs.append(os.path.join(os.environ['TRM_SOFTWARE'], 'lib'))
    include_dirs.append(os.path.join(os.environ['TRM_SOFTWARE'], 'include'))
else:
    print >>sys.stderr, "Environment variable TRM_SOFTWARE pointing to location of shareable libraries and includes not defined!"

include_dirs.append(numpy.get_include())

ext_modules = [
    Extension("lfit",
        ["lfit.pyx"],
        include_dirs = include_dirs,
        library_dirs = library_dirs,
        libraries = ["subs","roche"]
    )
]
    
setup(
    name = "lfit",
    ext_modules = cythonize(ext_modules)
)