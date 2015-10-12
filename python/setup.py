from __future__ import print_function
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy 
import os, sys

library_dirs = []
include_dirs = ['./include']

# need to direct to where includes and  libraries are
if 'TRM_SOFTWARE' in os.environ:
    library_dirs.append(os.path.join(os.environ['TRM_SOFTWARE'], 'lib'))
    include_dirs.append(os.path.join(os.environ['TRM_SOFTWARE'], 'include'))
else:
    print ("Environment variable TRM_SOFTWARE pointing to location of shareable libraries and includes not defined!",file=sys.stderr)
    sys.exit(-1)
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
    version="0.1",
    description="Calculate and fit CV lightcurves",
    ext_modules = cythonize(ext_modules),
    url = "https://github.com/StuartLittlefair/lfit",
    author_email = "s.littlefair@shef.ac.uk",
    scripts = ['scripts/calcPhysicalParams.py']
)