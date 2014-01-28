#include "trm/subs.h"
#include "trm/array1d.h"
#include "trm/array2d.h"
#include "trm/format.h"
#include "trm/input.h"
#include <cfloat>
#include <climits>
#include <ctime>
#include <iostream>
#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

using namespace std;

int
gsl_linalg_cholesky_decomp (gsl_matrix * A)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M != N)
    {
      GSL_ERROR("cholesky decomposition requires square matrix", GSL_ENOTSQR);
    }
  else
    {
      size_t i,j,k;
      int status = 0;

      /* Do the first 2 rows explicitly.  It is simple, and faster.  And
       * one can return if the matrix has only 1 or 2 rows.  
       */

      double A_00 = gsl_matrix_get (A, 0, 0);
      
      double L_00 = sqrt(A_00);
      
      if (A_00 <= 0)
        {
          status = GSL_EDOM ;
        }

      gsl_matrix_set (A, 0, 0, L_00);
  
      if (M > 1)
        {
          double A_10 = gsl_matrix_get (A, 1, 0);
          double A_11 = gsl_matrix_get (A, 1, 1);
          
          double L_10 = A_10 / L_00;
          double diag = A_11 - L_10 * L_10;
          double L_11 = sqrt(diag);
          
          if (diag <= 0)
            {
              status = GSL_EDOM;
            }

          gsl_matrix_set (A, 1, 0, L_10);        
          gsl_matrix_set (A, 1, 1, L_11);
        }
      
      for (k = 2; k < M; k++)
        {
          double A_kk = gsl_matrix_get (A, k, k);
          
          for (i = 0; i < k; i++)
            {
              double sum = 0;

              double A_ki = gsl_matrix_get (A, k, i);
              double A_ii = gsl_matrix_get (A, i, i);

              gsl_vector_view ci = gsl_matrix_row (A, i);
              gsl_vector_view ck = gsl_matrix_row (A, k);

              if (i > 0) {
                gsl_vector_view di = gsl_vector_subvector(&ci.vector, 0, i);
                gsl_vector_view dk = gsl_vector_subvector(&ck.vector, 0, i);
                
                gsl_blas_ddot (&di.vector, &dk.vector, &sum);
              }

              A_ki = (A_ki - sum) / A_ii;
              gsl_matrix_set (A, k, i, A_ki);
            } 

          {
            gsl_vector_view ck = gsl_matrix_row (A, k);
            gsl_vector_view dk = gsl_vector_subvector (&ck.vector, 0, k);
            
            double sum = gsl_blas_dnrm2 (&dk.vector);
            double diag = A_kk - sum * sum;

            double L_kk = sqrt(diag);
            
            if (diag <= 0)
              {
                status = GSL_EDOM;
              }
            
            gsl_matrix_set (A, k, k, L_kk);
          }
        }

      /* Now copy the transposed lower triangle to the upper triangle,
       * the diagonal is common.  
       */
      
      for (i = 1; i < M; i++)
        {
          for (j = 0; j < i; j++)
            {
              double A_ij = gsl_matrix_get (A, i, j);
              gsl_matrix_set (A, j, i, A_ij);
            }
        } 
      
      if (status == GSL_EDOM)
        {
	  std::string err = "matrix must be positive definite";
          //GSL_ERROR ("matrix must be positive definite", GSL_EDOM);
	  throw err;
        }
      
      return GSL_SUCCESS;
    }
}


int
gsl_linalg_cholesky_solve (const gsl_matrix * LLT,
                           const gsl_vector * b,
                           gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* Copy x <- b */

      gsl_vector_memcpy (x, b);

      /* Solve for c using forward-substitution, L c = b */

      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);

      /* Perform back-substitution, U x = c */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LLT, x);


      return GSL_SUCCESS;
    }
}

int
gsl_linalg_cholesky_svx (const gsl_matrix * LLT,
                         gsl_vector * x)
{
  if (LLT->size1 != LLT->size2)
    {
      GSL_ERROR ("cholesky matrix must be square", GSL_ENOTSQR);
    }
  else if (LLT->size2 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      /* Solve for c using forward-substitution, L c = b */

      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);

      /* Perform back-substitution, U x = c */

      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LLT, x);

      return GSL_SUCCESS;
    }
}

//template <class X>
Subs::Array2D<double> choleskyDecomp(const Subs::Array2D<double>& x){

  int nx = x.get_nx(), ny=x.get_ny();
  Subs::Array2D<double> y(nx,ny);

  gsl_matrix * m     = gsl_matrix_alloc (nx, ny);
  for(int i=0; i<x.get_ny(); i++){
    for(int j=0;j<x.get_nx();j++){
      gsl_matrix_set (m,     i, j, x[j][i]);
    }
  }

  try{
    gsl_linalg_cholesky_decomp (m);
  }catch(std::string err){
    throw err;
  }

  for(int i=0; i<x.get_ny(); i++){
    for(int j=0;j<x.get_nx();j++){
      if(j<=i){
	y[j][i] = gsl_matrix_get (m, i, j);
      }else{
	y[j][i] = 0.0;
      }
    }
  }
  return y; 
}
