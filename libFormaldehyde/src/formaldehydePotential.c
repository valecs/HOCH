#include "formaldehydePotential.h"
#include "formaldehydeProperties.h"
#include <gradient.h>
#include <assert.h>
#include <gsl/gsl_blas.h>

#define FORMALDEHYDE_GRADIENT_H 1e-4

/*TODO: eliminate asserts*/

extern void getpot(const int * natom, const double xmass[36], double xx[3][36], double * v, double * v1, double * v2, double * s);

//package the whole thing together
const potential potential_hoch_bare_nop = {
  .v = &potential_formaldehyde_nop,
  .g = &gradient_formaldehyde_nop
};
const potential * const hochPotential_nop = &potential_hoch_bare_nop;

const potential potential_hoch_bare = {
  .v = &potential_formaldehyde,
  .g = &gradient_formaldehyde
};
const potential * const hochPotential = &potential_hoch_bare;

static double HOCHgetCMDist(const gsl_vector *r, size_t offset);




/*
  returns the value of the formaldehyde potential at r while respecting the comment in
  h2co-mrci.f, which states:
    "note h1 and h2 are not equivalent, it is
    assumed that h2 is further away from the
    other atoms than h1"
  We do this by ensuring that h2 is indeed farther from the total center of mass than h1.

*/
double potential_formaldehyde(const gsl_vector * r){
  assert(r->size == 12);

  static const int natom = 4;
  static const double xmass[36]={
    1837.15   // H
    ,29156.95 // O
    ,21874.66 // C
    ,1837.15  // H
  }; //should initialize first 4 compnents and leave the others

  //double xx[3][36];
  static double v = 0, v1 = 0, v2 = 0, s = 0;
  static double xx[3][36]; //transpose of fortran's 36x3

  for (size_t coor = 0; coor < 3; coor++){
    for (size_t atom = 0; atom < 4; atom++){
      xx[coor][atom]=gsl_vector_get(r,(atom*3)+coor);
    }
  }

  if (HOCHgetCMDist(r,formaldehydeOffset.H1) > HOCHgetCMDist(r,formaldehydeOffset.H2)){
    double temp;
    const int H1 = 0;
    const int H2 = 3;
     
    for (size_t coor = 0; coor < 3; coor++){
      temp = xx[coor][H1];
      xx[coor][H1] = xx[coor][H2];
      xx[coor][H2] = temp;
    }
  }
  
  getpot(&natom, xmass, xx, &v, &v1, &v2, &s);

  return v;
}
/*
  gradient using the above function
*/
void gradient_formaldehyde(const gsl_vector * r, gsl_vector * grad){
  assert((r->size == 12) && (grad->size == 12));
  grad_5pt(&potential_formaldehyde, r, FORMALDEHYDE_GRADIENT_H, grad);
  return;
}

/*
  gives the distance of the offset-specified atom from the total center of mass
  use like:
  HOCHgetCMDist(r,formaldehydeOffset.H1) > HOCHgetCMDist(r,formaldehydeOffset.H2);
*/
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


/*
  returns the value of the formaldehyde potential at r without respecting permutation requirements
*/
double potential_formaldehyde_nop(const gsl_vector * r){
  assert(r->size == 12);

  static const int natom = 4;
  static const double xmass[36]={
    1837.15   // H
    ,29156.95 // O
    ,21874.66 // C
    ,1837.15  // H 
  }; //should initialize first 4 compnents and leave the others

  //double xx[3][36]; 
  static double v = 0, v1 = 0, v2 = 0, s = 0;
  static double xx[3][36]; //transpose of fortran's 36x3
  
  for (size_t coor = 0; coor < 3; coor++){
    for (size_t atom = 0; atom < 4; atom++){
      xx[coor][atom]=gsl_vector_get(r,(atom*3)+coor);
    }
  }
  
  getpot(&natom, xmass, xx, &v, &v1, &v2, &s);

  return v;
}

/*
  computes the gradient of the formaldehyde potential (no permutation) at r and stores the result in grad
 */
void gradient_formaldehyde_nop(const gsl_vector * r, gsl_vector * grad){
  assert((r->size == 12) && (grad->size == 12));
  grad_5pt(&potential_formaldehyde_nop, r, FORMALDEHYDE_GRADIENT_H, grad);
  return;
}
