#include "hessian.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

void compute_hessian(double (*v)(const gsl_vector *), gsl_matrix * H, const gsl_vector * x0,  double delta){
  size_t dim = x0 -> size;

  double f0, f1, f2, f3, f4;
  gsl_vector * ei, * ej, * x;
  
  ei = gsl_vector_alloc(dim);
  ej = gsl_vector_alloc(dim);
  x  = gsl_vector_alloc(dim);

  gsl_matrix_set_zero (H);

  // loop 1 takes care of the diagonal; see A&S 25.3.23
  f0 = (*v)(x0); // f0 is safe during this loop
  for (size_t i = 0; i < dim; i++){
    gsl_vector_set_basis(ei, i);

    gsl_vector_memcpy(x, x0);//dest,src
    gsl_blas_daxpy(delta,ei,x);// a,x,y ;y <- ax+y
    f1 = (*v)(x);

    gsl_blas_daxpy((-2.0)*delta,ei,x); //x0-2h
    f2 = (*v)(x);

    f3 = (1/(delta*delta))*(f1 + f2 - (2 * f0));

    gsl_matrix_set(H, i, i, f3);
  }

  // loop 2 does off-diagonal terms; see A&S 25.3.26
  // f0 destroyed in this loop
  for (size_t i = 0; i < dim; i++){
    gsl_vector_set_basis(ei, i);
    for (size_t j = 0 ; j < i; j++){
      gsl_vector_set_basis(ej, j);


      // don't thoughtlessly change the oder of opperations here; we go around
      // a square clockwise; see 2013.09.13:61

      gsl_vector_memcpy(x, x0);//dest,src

      gsl_blas_daxpy(delta, ei, x);
      gsl_blas_daxpy(delta, ej, x);
      f1 = (*v)(x);

      gsl_blas_daxpy((-2.0 * delta),ei,x);
      f2 = (*v)(x);

      gsl_blas_daxpy((-2.0 * delta),ej,x);
      f3 = (*v)(x);

      gsl_blas_daxpy((2.0 * delta),ei,x);
      f4 = (*v)(x);

      f0 = (1/(4*delta*delta))*(f1 + f3 - f2 - f4);

      gsl_matrix_set(H, i, j, f0);
      gsl_matrix_set(H, j, i, f0);

    }
  }

  gsl_vector_free(ei);
  gsl_vector_free(ej);
  gsl_vector_free(x);

}



#define M 5

/* Higher order hessian calculation, accurate to O(h^4) */
void hessian4(double (*F)(const gsl_vector *), gsl_matrix * H, const gsl_vector * x0, const double h){
  /* diagonal stencil */
  const int B[M] = {-1, 16, -30, 16, -1};
  /* off-diagonal stencil */
  const int C[M][M] = {{ 1,  -8, 0,  8,  -1},
		       {-8,  64, 0, -64,  8},
		       { 0,   0, 0,   0,  0},
		       { 8, -64, 0,  64, -8},
		       {-1,   8, 0,  -8,  1}};
  /* constant in denominator */
  const int D = 12;
  /* offset (since: central difference)  */
  const int p = 2;

  const size_t n = x0->size;
  gsl_matrix_set_zero(H);

  gsl_vector * x = gsl_vector_alloc(n);
  gsl_vector * ihat = gsl_vector_alloc(n);
  gsl_vector * jhat = gsl_vector_alloc(n);

  double * Hij = NULL;

  /* diagonal elements */
  for (size_t i = 0; i < n; i++){
    gsl_vector_set_basis(ihat, i);
    Hij = gsl_matrix_ptr(H, i, i);
    /* elements of stencil */
    for (int k = 0; k < M; k++){
      gsl_vector_memcpy(x, x0);
      gsl_blas_daxpy(h*(k-p), ihat, x);
      (*Hij) += (*F)(x) * B[k];
    }
    (*Hij) /= D * h * h;
  }

  /* lower off-diagonal elements */
  for (size_t i = 0; i < n; i++){
    for (size_t j = 0; j < i; j++){
      gsl_vector_set_basis(ihat, i);
      gsl_vector_set_basis(jhat, j);
      Hij = gsl_matrix_ptr(H, i, j);
      /* elements of stencil */
      for (int k = 0; k < M; k++){
	for (int l = 0; l < M; l++){
	  if (!C[k][l]){continue;} // skip 0 elements
	  
	  gsl_vector_memcpy(x, x0);
	  gsl_blas_daxpy(h*(k-p), ihat, x);
	  gsl_blas_daxpy(h*(l-p), jhat, x);
	  (*Hij) += (*F)(x) * C[k][l];
	}
      }
      (*Hij) /= D * D * h * h;
      /* set corresponding upper element */
      gsl_matrix_set(H, j, i, (*Hij));
    }
  }
  
  gsl_vector_free(x);
  gsl_vector_free(ihat);
  gsl_vector_free(jhat);
}
