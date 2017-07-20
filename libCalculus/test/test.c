#include "../src/gradient.h"
#include "../src/hessian.h"

#include <potentials.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <math.h>

static const size_t dim = 2;

static double bowl(const gsl_vector * r);
static int eq(double a, double b, double tolerance);

int main (void){
  energy v = &bowl;
  gsl_vector * grad = gsl_vector_alloc(dim);
  gsl_vector * r = gsl_vector_alloc(dim);
  gsl_matrix * H = gsl_matrix_alloc(dim,dim);

  const double x1 = 2.0;
  const double x2 = 3.0;
  const double tol = 1e-6;
  const double delta = 1e-4;

  gsl_vector_set(r, (size_t) 0, x1);
  gsl_vector_set(r, (size_t) 1, x2);

  if (GSL_SUCCESS != grad_5pt(v, r, delta, grad)){
    return EXIT_FAILURE;
  }

  if (!eq(gsl_vector_get(grad, (size_t) 0), 2*x1, tol) ||
      !eq(gsl_vector_get(grad, (size_t) 1), 2*x2, tol)){
  return EXIT_FAILURE;
  }

  compute_hessian(v, H, r, delta);

  if (!eq(gsl_matrix_get(H, (size_t) 0, (size_t) 0), 2.0, tol) ||
      !eq(gsl_matrix_get(H, (size_t) 0, (size_t) 1), 0.0, tol) ||
      !eq(gsl_matrix_get(H, (size_t) 1, (size_t) 0), 0.0, tol) ||
      !eq(gsl_matrix_get(H, (size_t) 1, (size_t) 1), 2.0, tol)){
  return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

static double bowl(const gsl_vector * r){
  double v = 0;
  gsl_blas_ddot(r,r, &v);
  return v;
}

static int eq(double a, double b, double tolerance){
  if (fabs(a-b) < tolerance)
    return 1;
  else
    return 0;
}
