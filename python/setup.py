from __future__ import print_function
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy
import os, sys
import platform

library_dirs = []
include_dirs = ["./include"]

# need to direct to where includes and  libraries are
if "TRM_SOFTWARE" in os.environ:
    library_dirs.append(os.path.join(os.environ["TRM_SOFTWARE"], "lib"))
    include_dirs.append(os.path.join(os.environ["TRM_SOFTWARE"], "include"))
else:
    print(
        "Environment variable TRM_SOFTWARE pointing to location of shareable libraries and includes not defined!",
        file=sys.stderr,
    )
    sys.exit(-1)
include_dirs.append(numpy.get_include())

if platform.system() == "Darwin":
    ext_modules = [
        Extension(
            "lfit",
            ["lfit.pyx"],
            include_dirs=include_dirs,
            library_dirs=library_dirs,
            extra_compile_args=["-stdlib=libc++"],
            extra_link_args=["-stdlib=libc++"],
            language="c++",
            libraries=["subs", "roche"],
        )
    ]
else:
    ext_modules = [
        Extension(
            "lfit",
            ["lfit.pyx"],
            include_dirs=include_dirs,
            library_dirs=library_dirs,
            language="c++",
            libraries=["subs", "roche"],
        )
    ]

setup(
    name="lfit",
    version="0.15",
    description="Calculate and fit CV lightcurves",
    ext_modules=cythonize(ext_modules),
    url="https://github.com/StuartLittlefair/lfit",
    author_email="s.littlefair@shef.ac.uk",
)
