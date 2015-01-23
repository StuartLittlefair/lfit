README ======

LFIT is a C++ program for fitting the lightcurves of eclipsing Dwarf
Nova. It leans *very* heavily on Tom Marsh's C++ library for Roche
geometry and it's dependencies.

LFIT also have a python (cython) wrapper to allow the computation and
fitting of CV lightcurves in the Python programming language.

INSTALLATION ------------

#C++ program
To install the C++ program, edit the Makefile until the line
STOP_EDITING. Then simply type

 make lfit

to compile the complete C++ program, or

 make lfit_fake

to make a simple program that just computes lightcurves

#params
LFIT uses a set of 14 parameters to describe a CV. These can be combined
with White-Dwarf mass radius measurements to derive the mass and radius
of both components in the CV. Routines for this are in the params subdirectory.

#PYTHON BINDINGS
The python binding is in the sub-directory "python". As well as numpy,
you will also need cython installed to create the wrapper.

Installation proceeds via the usual::

 python setup.py install

if you are root, or::

 python setup.py install --prefix=<install dir>

if you are not.

If you want to test the module use

 python setup.py build_ext --inplace


KNOWN_ISSUES -------------

None yet.



