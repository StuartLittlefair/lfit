#include "trm/subs.h"
#include "trm/array2d.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

int gsl_linalg_cholesky_decomp(gsl_matrix * A);

//template <class X>
Subs::Array2D<double> choleskyDecomp(const Subs::Array2D<double>& x);

