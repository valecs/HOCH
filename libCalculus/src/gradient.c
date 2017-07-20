#include "gradient.h"

#include "potentials.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

//returns an approximation, via a 5-point stencil with spacing h along each axis, of the gradient of f at x0
int grad_5pt(energy f, const gsl_vector * x0, double h, gsl_vector * grad){
  size_t dim = x0 -> size;

  static gsl_vector * xihat = NULL;
  static gsl_vector * point = NULL;
  
  if (xihat == NULL){xihat = gsl_vector_alloc(dim);}
  if (point == NULL){point = gsl_vector_alloc(dim);}

  if ((xihat->size != dim) ||
      (point->size != dim)){
    return GSL_EBADLEN;
  }

  //if doing the following manually, would have something like:
  /*
    volitile double temp; 
    temp = xoi+h;
    h = temp-xio;
  */
  //on each loop


  double dfdxi;
  
  for (size_t i = 0; i < dim; i++){
  
    gsl_vector_set_basis(xihat, i);
    dfdxi=0;

    //+f(x - 2h)/12
    gsl_vector_memcpy(point, x0);
    gsl_blas_daxpy((-2.0)*h,xihat,point);
    dfdxi += ((*f)(point))/12.0;

    //-f(x - h) * 2/3
    gsl_vector_memcpy(point, x0);
    gsl_blas_daxpy((-1.0)*h,xihat,point);
    dfdxi -= 2.0*((*f)(point))/3.0;

    //-f(x + h) * 2/3
    gsl_vector_memcpy(point, x0);
    gsl_blas_daxpy((1.0)*h,xihat,point);
    dfdxi += 2.0*((*f)(point))/3.0;

    //-f(x + 2h)/12
    gsl_vector_memcpy(point, x0);
    gsl_blas_daxpy((2.0)*h,xihat,point);
    dfdxi -= ((*f)(point))/12.0;

    gsl_vector_set(grad, i, dfdxi);
  }

  gsl_vector_scale(grad, 1/h);

  //gsl_vector_free(xihat); //since we carry them around now
  //gsl_vector_free(point); //since we carry them around now

  return GSL_SUCCESS;
}
