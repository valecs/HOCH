#include "../src/formaldehydePotential.h"
#include "../src/formaldehydeProperties.h"
#include "../src/formaldehydeRoaming.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

double HOCHgetCMDist(const gsl_vector *r, size_t offset);

static const double geometry_array[12]={
  -1.12460784, 1.79389857, 0.0
  , 2.2882014, 0.0, 0.0
  , 0.0, 0.0, 0.0
  , -1.12460784, -1.79389857, 0.0
};

static const size_t dim = 12;

void roamTest(void);

int main(void){
  gsl_vector_const_view geometry_view = gsl_vector_const_view_array (geometry_array, dim);
  double value = (*(hochPotential->v))(&geometry_view.vector);
  printf("Got %#15.15g\n", value);
  
  //roamTest();
  /*
    run with the above and:
    ./test < ~/src/formaldehyde/analysis/geodesics/roamingRegion/5a0/region-5a0.points \
      | wc \
      | awk '{ print $1"/12" }' \
      | bc -l

    should yield 1872
  */

  //cmDistTest
  for (int i = 0; i < formaldehydeNumCenters; i++){
    printf("%#15.15g\n",HOCHgetCMDist(&geometry_view.vector, formaldehydeOffsetArray[i]));
  }
  
  if (value < 1e-10)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

void roamTest(void){
  gsl_vector * v = gsl_vector_alloc(dim);

  gsl_error_handler_t * err_gsl = gsl_set_error_handler_off(); //restore with gsl_set_error_handler(err_gsl);

  while(gsl_vector_fscanf(stdin, v) == GSL_SUCCESS){
    if (isRoaming(v))
      gsl_vector_fprintf(stdout, v, "%#15.15g");
  }
  gsl_set_error_handler(err_gsl);

  return;
}

double HOCHgetCMDist(const gsl_vector *r, size_t offset){
  const size_t space = 3;

  double cm_array[space];
  gsl_vector_view cm = gsl_vector_view_array(cm_array, space);
  gsl_vector_set_zero(&cm.vector);
  
  for (int i = 0; i < formaldehydeNumCenters; i++){
    gsl_vector_const_view v = gsl_vector_const_subvector(r,formaldehydeOffsetArray[i],space);
    gsl_blas_daxpy(
		   formaldehydeMassArray[i] / formaldehydeMassTotal
		   , &v.vector
		   , &cm.vector);
  }

  gsl_vector_const_view atom = gsl_vector_const_subvector(r, offset, space);
  gsl_vector_sub(&cm.vector, &atom.vector);
  
  return gsl_blas_dnrm2(&cm.vector);
}
