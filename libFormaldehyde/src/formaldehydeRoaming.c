#include "formaldehydeRoaming.h"
#include "formaldehydePotential.h"
#include "formaldehydeProperties.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

static void HOCHgetCM(const gsl_vector *r, gsl_vector * cm);
static double dist(const gsl_vector * a, const gsl_vector * b);

/*
Returns:
   1 , if the configuration space vector r satisfies the roaming criteria of 2014.08.14:174
   0 , othwise

Briefly, the criteria are:
  [Min Force]: F(H-CM) > -2.5x10^-4 (EH/a0)
  [Min Distance]: |rHH| > 5.86728 a0
  [Max Distance]: |rHH| < 9.0 a0
*/

int isRoaming(const gsl_vector * r){
  const double minForce = -2.5e-4;
  const double minDistance = 5.86728;
  const double maxDistance = 9.0;

  const size_t space = 3;
  const struct offset_t o = formaldehydeOffset;

  const int presentlyRoaming = 1, notRoaming = 0;

  /*TODO: free these*/
  static gsl_vector * grad = NULL; //gradient
  static gsl_vector * cm = NULL; // center of mass
  static gsl_vector * rH = NULL; // H position
  //static gsl_vector * fH = NULL; // H force

  if (!grad){grad =  gsl_vector_alloc(formaldehydeNumDimensions);}
  if (!cm){cm = gsl_vector_alloc(space);}
  if (!rH){rH = gsl_vector_alloc(space);}
  //  if (!fH){fH = gsl_vector_alloc(space);}

  gsl_vector_const_view rH1 =  gsl_vector_const_subvector(r, o.H1, space);
  gsl_vector_const_view rH2 =  gsl_vector_const_subvector(r, o.H2, space);

  double distance = dist(&rH1.vector,&rH2.vector);

  if ((distance < minDistance) ||
      (distance > maxDistance)){
    return notRoaming;
  }


  size_t roamerOffset;
  HOCHgetCM(r,cm);

  if (dist(cm,&rH1.vector) > dist(cm,&rH2.vector))
    roamerOffset = o.H1;
  else
    roamerOffset = o.H2;

  //consider factoring out the below for increased clarity/functionality at the expense of some repeat computation

  gsl_vector_const_view rHr = gsl_vector_const_subvector(r, roamerOffset, space);

  gsl_vector_memcpy(rH, &rHr.vector);
  gsl_vector_sub(rH, cm);
  gsl_vector_scale(rH, -1.0/(gsl_blas_dnrm2(rH))); // yeilds the negative unit vector since we don't take the negative of the gradient.

  (*(hochPotential->g))(r,grad);
  gsl_vector_view fH = gsl_vector_subvector(grad, roamerOffset, space);

  double forceOnRoamingH = 0;
  gsl_blas_ddot(rH,&fH.vector, &forceOnRoamingH);

  if (forceOnRoamingH < minForce){
    return notRoaming;
  }

  // making it here means we roam!

  return presentlyRoaming;
}

// specialized function for formaldehyde cm
static void HOCHgetCM(const gsl_vector *r, gsl_vector * cm){
  const size_t space = 3;
  gsl_vector_set_zero(cm);
  for (int i = 0; i < formaldehydeNumCenters; i++){
    gsl_vector_const_view v = gsl_vector_const_subvector(r,formaldehydeOffsetArray[i],space);
    gsl_blas_daxpy(
		   formaldehydeMassArray[i] / formaldehydeMassTotal
		   , &v.vector
		   , cm);
  }
  return;
}

//returns the distance between 2 vectors
static double dist(const gsl_vector * a, const gsl_vector * b){
  gsl_vector * t = gsl_vector_alloc(a->size);
  double d = 0;
  gsl_vector_memcpy(t,a);

  gsl_vector_sub(t,b);
  d = gsl_blas_dnrm2(t);

  gsl_vector_free(t);

  return d;
}
