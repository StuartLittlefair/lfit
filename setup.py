from __future__ import print_function

import os
import platform
import sys

import numpy
from Cython.Build import cythonize
from setuptools import Extension, setup

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
    # MacOS specific settings
    extra_compile_args = ["-stdlib=libc++"]
    extra_link_args = ["-stdlib=libc++"]
else:
    extra_compile_args = []
    extra_link_args = []


extensions = [
    Extension(
        "lfit",
        [
            "lfit.pyx",
            "src/WhiteDwarf.cc",
            "src/Disc.cc",
            "src/BrightSpot.cc",
            "src/Donor.cc",
            "src/finddeg.cc",
            "src/Point.cc",
        ],
        include_dirs=include_dirs,
        library_dirs=library_dirs,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        language="c++",
        libraries=["subs", "roche"],
    )
]

setup(
    name="lfit",
    version="0.15",
    description="Calculate and fit CV lightcurves",
    ext_modules=cythonize(extensions, language_level="3"),
    url="https://github.com/StuartLittlefair/lfit",
    author_email="s.littlefair@shef.ac.uk",
)
