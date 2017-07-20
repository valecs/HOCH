#include "HOCHgeodesics.h"

#include <formaldehydeProperties.h>
#include <formaldehydePotential.h>
#include <linkedList.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>

#include <math.h> //for: double ceil(double x);
#include <stdbool.h>

//static (gsl_vector *)(*dup)(gsl_vector *) = &gsl_vector_duplicate;

static llist * shortestDestroyLong(llist * A, llist * B);
static int burrowWithMasses(gsl_vector * r, double landscapeEnergy);
static double HOCHkinematicPathLength(llist * p);
static llist * generatePerturbedPath(llist * path, double landscapeEnergy, gsl_rng * r);
static gsl_vector * gsl_vector_duplicate(const gsl_vector * v);
static double gsl_vector_dist(const gsl_vector * u, const gsl_vector * v);
static void getCM(const gsl_vector * v, gsl_vector * cm);

/* declared in header:

llist * getOptimizedHOCHGeodesic(const gsl_vector * A, const gsl_vector * B, double landscapeEnergy);
llist * findHOCHGeodesic(const gsl_vector * A, const gsl_vector * B, double landscapeEnergy);
llist * optimizeHOCHGeodesic(llist * path, double landscapeEnergy);
void path_free(llist * path);
void resetCenterOfMass(gsl_vector * v);

 */

llist * getOptimizedHOCHGeodesic(const gsl_vector * A, const gsl_vector * B, double landscapeEnergy){
  llist * ABpath = findHOCHGeodesic(A, B, landscapeEnergy);
  llist * BApath = findHOCHGeodesic(B, A, landscapeEnergy);
  BApath = reverse(BApath); //make both lists have the same order

  ABpath = optimizeHOCHGeodesic(ABpath, landscapeEnergy);
  BApath = optimizeHOCHGeodesic(BApath, landscapeEnergy);

  return shortestDestroyLong(ABpath, BApath);
}

/*Returns the path of shortest kinematic length
  and frees the longer one.
*/
static llist * shortestDestroyLong(llist * A, llist * B){
  if (!A && !B)
    return NULL;
  else if (!A)
    return B;
  else if (!B)
    return A;

  if (HOCHkinematicPathLength(A) < HOCHkinematicPathLength(B)){
    path_free(B);
    return A;
  } else{
    path_free(A);
    return B;
  }
}

llist * findHOCHGeodesic(const gsl_vector * A, const gsl_vector * B, double landscapeEnergy){
  gsl_vector * origin  = gsl_vector_duplicate(A);
  gsl_vector * target = gsl_vector_duplicate(B);

  llist * path = ll_new(origin);

  energy v = hochPotential->v;
  const double deltaS = 5e-2;
  
  // changed the scaling constant 1000->100, see 2015.01.08:14
  const int maxItr = (int) ceil((gsl_vector_dist(origin,target)/deltaS) * 100);

  gsl_vector * r = gsl_vector_duplicate(origin);
  gsl_vector * dr = gsl_vector_alloc(formaldehydeNumDimensions);

  int itr = 0;
  while(gsl_vector_dist(r, target) > 2 * deltaS){
    int burrowStatus = GSL_SUCCESS;

    gsl_vector_memcpy(dr,target);
    gsl_vector_sub(dr,r); //dr <- rf - rt
    double magDr = gsl_blas_dnrm2(dr);

    //step along
    gsl_blas_daxpy(deltaS/magDr,dr,r);

    //burrow if necessary
    if ((*v)(r) > landscapeEnergy)
      burrowStatus = burrowWithMasses(r, landscapeEnergy);

    if ((burrowStatus != GSL_SUCCESS) ||
	(itr++ > maxItr)){
	path_free(path); path = NULL;
	break;
    }

    path = append(path, gsl_vector_duplicate(r));
  }
  if(path) //only add the last point if we didn't fail
    path = append(path, target);
  else
    gsl_vector_free(target);

  gsl_vector_free(r);
  gsl_vector_free(dr);

  return path; //can be NULL!
}


/* Computes a burrowing trajectory using the center of mass invariant projection described in
   "Center of Mass Translation in the Geodesic", 2014.10.27
*/
static int burrowWithMasses(gsl_vector * r, double landscapeEnergy){
  energy v = hochPotential->v;
  gradient g = hochPotential->g;

  gsl_matrix_const_view CMIP = gsl_matrix_const_view_array (formaldehydeCMInvariantProjectionArray
							    ,formaldehydeNumDimensions
							    ,formaldehydeNumDimensions);
  int status = GSL_SUCCESS;
  const int maxItr = 500;
  const double fractionalTolerance = 1e-6;

  gsl_vector * grad = gsl_vector_alloc(formaldehydeNumDimensions);

  /* 
     Guarantee at least one iteration since we enter this function from a V(R) > E_L condition
  */

  int itr = 0;
  do{
    if (itr++ > maxItr){
      status = GSL_FAILURE;
      break;
    }

    /*
      int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x
                        , double beta, gsl_vector * y);

      Compute the matrix-vector product and sum:
      y = \alpha op(A) x + \beta y,
      where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans. 

      Here:
      r <- (alpha * CMIP * grad) + r
    */

    (*g)(r,grad); //gradient
    
    double gradMag2;
    gsl_blas_ddot(grad, grad, &gradMag2);
    double alpha = -((*v)(r) - landscapeEnergy)/gradMag2;

    status = gsl_blas_dgemv((CBLAS_TRANSPOSE_t) CblasNoTrans, alpha, &CMIP.matrix, grad, 1.0, r);

    /*{//debugging block
      double temp_array[formaldehydeNumDimensions];
      gsl_vector_view temp = gsl_vector_view_array(temp_array, formaldehydeNumDimensions);

      double temp3_array[3];
      gsl_vector_view temp3 = gsl_vector_view_array(temp3_array, (size_t) 3);

      double m_corrected, m_old, m_spurrious;

      gsl_blas_dgemv((CBLAS_TRANSPOSE_t) CblasNoTrans, alpha, &CMIP.matrix, grad, 0.0, &temp.vector);
      getCM(&temp.vector, &temp3.vector);
      m_corrected = gsl_blas_dnrm2(&temp3.vector);

      getCM(grad, &temp3.vector);
      m_old = -1.0 * alpha * gsl_blas_dnrm2(&temp3.vector);
      m_spurrious = 0.0;

      //fprintf(stdout,"%17.14g\t%17.14g\t%17.14g\n", m_corrected, m_old, m_spurrious);
      fprintf(stdout,"%17.14g + ", m_corrected);
      }*/
  }while((*v)(r) > landscapeEnergy * (1+fractionalTolerance));
  
  gsl_vector_free(grad);

  /*
    debugging block, goal:
      meausre the shift in center of mass, i.e. before and after the reset to see if the
      many v. small, O(1e-19), deviations are compounding crazilly, somehow.
  */  

  /*  
  double pre_array[3];
  double post_array[3];
  gsl_vector_view pre  = gsl_vector_view_array(pre_array, (size_t) 3);
  gsl_vector_view post = gsl_vector_view_array(post_array, (size_t) 3);

  getCM(r,&pre.vector);
  */

  //This is a shim to see if we manage to get the proper behavior out; I'll go back and fix the math if this works. If it does not, it isolates the problem away from here.
  //resetCenterOfMass(r);

  /*
  getCM(r,&post.vector);

  double xlation = gsl_vector_dist(&pre.vector, &post.vector);
  fprintf(stdout, "-> %17.14g\n", xlation);
  */
  return status;
}

static double HOCHkinematicPathLength(llist * p){
  if (!p)
    return NAN;//empty path has undefined length;
  
  gsl_matrix_const_view M = gsl_matrix_const_view_array(formaldehydeMassMatrixArray, formaldehydeNumDimensions, formaldehydeNumDimensions);

  const size_t dim = formaldehydeNumDimensions;

  double sum = 0, res;

  gsl_vector * rdot = gsl_vector_alloc(dim);
  gsl_vector * temp = gsl_vector_alloc(dim);

  p = first(p);

  //for(p=first(p); (p && p->succ); p=p->succ){
  while(p->succ != NULL){

    gsl_vector_memcpy(temp, (gsl_vector*) p->e); // temp <- xi
    gsl_vector_memcpy(rdot, (gsl_vector*) p->succ->e); // rdot <- xi+1

    gsl_vector_sub(rdot,temp); //rdot <- xi+1 - xi = rdot
    gsl_vector_memcpy(temp,rdot); //temp <- xi+1 - xi = rdot

    //rdot <- M . (xi+1 - xi)
    gsl_blas_dgemv ((CBLAS_TRANSPOSE_t) CblasNoTrans
		    , 1.0, &M.matrix, rdot, 0.0, temp);

    gsl_blas_ddot(rdot,temp,&res); // res <- (xi+1 - xi) . M . (xi+1 - xi)

    sum += sqrt(res);

    p = p->succ;
  }

  gsl_vector_free(rdot);
  gsl_vector_free(temp);

  return sum;
}

/*
Steps:
  do
    pick perturbation element of path
    perturb
    compute geodesics with new boundary conditions
    if they yeild a shorter path, return the new one
  while convergence not met
*/
llist * optimizeHOCHGeodesic(llist * path, double landscapeEnergy){
  if (!path)
    return path; //nothing to optimize

  /*
    Keep track of which round we're on.
    This can be used to id the path being optimized

    e.g.: 1 & 2 refer to the forward and reverse optimzation of the 1st path
   */
  static int optRound = 0;
  optRound++;

  const int convergedItr = 1000;
  int itrSinceChange = 0;
  double lengthCurrent, lengthPrime;

  lengthCurrent = HOCHkinematicPathLength(path);

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937); //freed at eof

  do{
    llist * pathPrime = generatePerturbedPath(path, landscapeEnergy, r); //may return NULL
    lengthPrime = HOCHkinematicPathLength(pathPrime);                    //may, therefore, return NaN

    if (pathPrime && lengthPrime < lengthCurrent){             //must, therefore, compare in ths way 
      path_free(path); path = NULL;
      path = pathPrime;
      lengthCurrent = lengthPrime;
      //fprintf(stderr, "optimizing%3d\t%5d\t%16.16g\n", optRound, itrSinceChange, lengthCurrent);
      itrSinceChange = 0;
    }
    else{
      path_free(pathPrime); pathPrime = NULL;
      itrSinceChange++;
    }
  }while(itrSinceChange < convergedItr);

  gsl_rng_free (r);
  return path;
}

/*
  Takes in a path and returns a perturbed one:
   At a random point in the path, 
   some center is displaced by a random amount along some axis
   And new geodesics are computed
*/
static llist * generatePerturbedPath(llist * path, double landscapeEnergy, gsl_rng * r){
  int pathLength = (int) length(path);
  energy v = hochPotential->v;

  const double deltaP = 5e-2;

  //this guarntee is important to prevent over-wrappping when casting (pathLength - 2) to unsigned
  if (pathLength < 3) //need at least 1 element not on the boundary
    return NULL;

  gsl_vector * educt   = (gsl_vector *)(first(path)->e);
  gsl_vector * product = (gsl_vector *)(last (path)->e);
  gsl_vector * perturbation    = gsl_vector_alloc(formaldehydeNumDimensions);
  gsl_vector * randomDirection = gsl_vector_alloc(formaldehydeNumDimensions);

  gsl_matrix_const_view CMIP = gsl_matrix_const_view_array (formaldehydeCMInvariantProjectionArray
							    ,formaldehydeNumDimensions
							    ,formaldehydeNumDimensions);

  /*
    Need to select a random step in the path, which is not an endpoint.
    unsigned long int gsl_rng_uniform_int (const gsl_rng * r, unsigned long int n)
      returns a random integer on [0,n)
    llist * index(llist * curr, unsigned int ind)
      begins counting at 0

    We therefore need to sample a range with l-2 elements which begins with 1.
    i.e.: [1,l-2], which we can get from:
    1 + gsl_rng_uniform_int(r,l-2)
  */
  do{
    unsigned int randomIndex = 1 + (unsigned int) gsl_rng_uniform_int (r, (unsigned long int) pathLength - 2);
    gsl_vector_memcpy(perturbation, (gsl_vector *)(lindex(path,randomIndex)->e));


    /*
      FIXME: we should use the standard and simpler to implement x' = x + e*(d*4*(X-0.5))
      where d is \delta_p, X is a uniform random variable on [0,1] and e is a unit vector along a random coordinate
    */
    
    gsl_vector_set_basis(randomDirection,
			 (size_t) gsl_rng_uniform_int (r, formaldehydeNumDimensions)
			 );
    int randomSign = (gsl_rng_uniform(r) > 0.5) ? 1 : -1;

    //perturbation <- perturbation +/- dp*(1-M)*e
    gsl_blas_dgemv((CBLAS_TRANSPOSE_t) CblasNoTrans, deltaP * randomSign, &CMIP.matrix, randomDirection, 1.0, perturbation);
    //gsl_blas_daxpy(deltaP * randomSign, randomDirection, perturbation); //fails to preserve CM invariance

  }while((*v)(perturbation) > landscapeEnergy);
  gsl_vector_free(randomDirection);

  // FIXME: Should compute geodesic in both directions and take the shorter one
  llist *pathA, *pathB;
  pathA = findHOCHGeodesic(educt,   perturbation, landscapeEnergy);

  if (pathA)
    pathB = findHOCHGeodesic(perturbation, product, landscapeEnergy);
  else
    return NULL;

  if (pathA && pathB){
    gsl_vector_free(perturbation);
    perturbation = NULL;
    gsl_vector_free(first(pathB)->e);
    first(pathB)->e = NULL;
    return join(pathA, tail(pathB));
  }
  else{
    path_free(pathA);
    path_free(pathB);
    return NULL;
  }

}

static gsl_vector * gsl_vector_duplicate(const gsl_vector * v){
  gsl_vector * u = gsl_vector_alloc(v->size);
  gsl_vector_memcpy(u,v);
  return u;
}

//sets the center of mass of v to that used in this body of work.
void resetCenterOfMass(gsl_vector * v){
  const size_t space = 3;

  //target CM hardcoded from 2014.09.30:187
  static const double rT_array[3] = {1.1440231452,0,0};
  gsl_vector_const_view rT =  gsl_vector_const_view_array (rT_array, space);

  static double cm_array[3];
  gsl_vector_view cm = gsl_vector_view_array (cm_array, space);

  getCM(v, &cm.vector);
  gsl_vector_sub(&cm.vector, &rT.vector);

  for (int i = 0; i < formaldehydeNumCenters; i++ ){
    gsl_vector_view center = gsl_vector_subvector(v, formaldehydeOffsetArray[i], space);
    gsl_vector_sub(&center.vector, &cm.vector);
  }
}

//returns the center of mass of the vector v in cm
static void getCM(const gsl_vector * v, gsl_vector * cm){
  const size_t space = 3;
  gsl_vector_set_zero(cm);

  for (int i = 0; i < formaldehydeNumCenters; i++ ){
    gsl_vector_const_view center = gsl_vector_const_subvector(v, formaldehydeOffsetArray[i],space);
    gsl_blas_daxpy(formaldehydeMassArray[i]/formaldehydeMassTotal, &center.vector, cm);
  }
}

//Returns the euclidian distance between u and v
static double gsl_vector_dist(const gsl_vector * u, const gsl_vector * v){
  gsl_vector * r = gsl_vector_duplicate(u);
  gsl_vector_sub(r,v);
  double d = gsl_blas_dnrm2(r);
  gsl_vector_free(r);
  return d;  
}

void path_free(llist * path){
  if (!path)
    return; //null paths need not be freed

  path = first(path);
  while (path->succ != NULL){
    gsl_vector_free(path->e); path->e = NULL;
    path = path->succ;        
    free(path->pred);         path->pred = NULL;
  }
  gsl_vector_free(path->e);   path->e = NULL;
  free(path);
}
