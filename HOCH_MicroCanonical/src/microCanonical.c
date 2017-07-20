/*
  This file enables the generation of microcannonically sampled initial-
  conditions on potential surfaces when starting from a well. The
  procedure followed is taken from Hase in Advances in Chemical Physics,
  105. And uses a normal mode expansion

  Notes expanding on the work are found in 2013.09.11:52-2013.09.12:60
 */

#include "microCanonical.h"
#include <hessian.h>
#include <potentials.h>
#include "cmHelpers.h"
#include "pathHelpers.h"

#include "debug_print.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>

#include <math.h>
#include <assert.h>

//not defined in math.h for c99 mode
#define PI 3.14159265358979323846264338327

void construct_inverse_root_mass_matrix(gsl_vector * mass, gsl_matrix * rootmassinv);

int get_phase_space_point(gsl_vector  * x0   // initial position
			   ,gsl_vector * mass // vector of masses
			   ,energy       v    // potential
			   ,double       E0   // total desired energy
			   ,gsl_vector * phase_space_point
			   ,double       hessian_tolerance
			   ,double       energy_tolerance
			   ,double       angular_momentum_tolerance
			   ,double       cm_velocity_tolerance){
  size_t dim, modes;
  dim = x0 -> size;

  // make sure the dimensions are ok
  assert((dim == mass -> size) && (2 * dim == phase_space_point-> size));
  // make sure the potential gives reasonable values at x0
  assert(((*v)(x0) != (1.0/0.0)) && ((*v)(x0) != (-1.0/0.0)) && ((*v)(x0) != (0.0/0.0)));

  // allocate instance of "Mersenne Twister" PRNG
  // RANDLUX and TAUS algorithms are attractive as well; See GSL docs for more info.
  static gsl_rng * r = NULL;
  if (r == NULL){
    r = gsl_rng_alloc(gsl_rng_mt19937);
    /*
      never freed! We do this so we maintain the state of the RNG
      across calls. The alternative was passing the rng in each time;
      more sensicle would be this + passing a seed. Seed is set in
      main via: extern unsigned long int gsl_rng_default_seed;
    */
    
    //"warm up" sequence
    for (int i = 0; i<1e6; i++){gsl_rng_uniform_pos(r);}
  }

  gsl_vector * f; // freed at end
  f = gsl_vector_calloc(dim); // frequencies

  gsl_matrix * H, * M2v, * L, * temp; // H freed after use, others at end
  H    = gsl_matrix_calloc(dim, dim); // Hessian
  M2v  = gsl_matrix_calloc(dim, dim); // M^(-1/2)
  L    = gsl_matrix_calloc(dim, dim); // Normal coordinate transform matrix
  temp = gsl_matrix_calloc(dim, dim); // temp

  // Allocate GSL eigen workspace
  gsl_eigen_symmv_workspace * e;
  e = gsl_eigen_symmv_alloc(dim); //freed after use


  // Get M^(-1/2)
  construct_inverse_root_mass_matrix(mass, M2v);

  /* printf("M2v:\n"); */
  /* print_matrix(M2v); */

  // Generate the Hessian & weight by mass
  compute_hessian(v, H, x0, hessian_tolerance);

  // weighting be mass is a 3-stepper
  // H <- M^(1/2) H M^(1/2)
  gsl_matrix_set_zero(temp);
  gsl_matrix_memcpy(temp, H); //dest, src
  gsl_blas_dsymm ((CBLAS_SIDE_t) CblasLeft,  (CBLAS_UPLO_t) CblasUpper, 1.0, M2v, temp , 0.0, H);

  gsl_matrix_memcpy(temp, H);
  gsl_blas_dsymm ((CBLAS_SIDE_t) CblasRight, (CBLAS_UPLO_t) CblasUpper, 1.0, M2v, temp , 0.0, H);

  // Diagonalize & sort
  assert(GSL_SUCCESS == gsl_eigen_symmv(H, f, L, e));
  assert(GSL_SUCCESS == gsl_eigen_symmv_sort(f, L, (gsl_eigen_sort_t) GSL_EIGEN_SORT_VAL_DESC));


  /* //some tests */

  /* //construct L^THL */
  /* //  gsl_blas_dgemm ((CBLAS_TRANSPOSE_t) CblasNoTrans, (CBLAS_TRANSPOSE_t) CblasNoTrans , 1.0, L,    H, 0.0, temp); */
  /* gsl_blas_dgemm((CBLAS_TRANSPOSE_t) CblasTrans */
  /* 		 , (CBLAS_TRANSPOSE_t) CblasNoTrans */
  /* 		 , 1.0 */
  /* 		 , L */
  /* 		 , H */
  /* 		 , 0.0 */
  /* 		 , temp); */


  /* printf("LH:\n"); */
  /* print_matrix(temp); */

  /* gsl_blas_dgemm ((CBLAS_TRANSPOSE_t) CblasNoTrans, (CBLAS_TRANSPOSE_t) CblasNoTrans   , 1.0, temp, L, 0.0, H); */

  /* printf("L^THL:\n"); */
  /* print_matrix(H); */

  /* //consrtuct LL^T, store result in H, since it's about to be freeed  */
  /* gsl_matrix_memcpy(temp, L); //dest, src   */
  /* gsl_blas_dgemm ((CBLAS_TRANSPOSE_t) CblasNoTrans, (CBLAS_TRANSPOSE_t) CblasTrans , 1.0, L, temp, 0.0, H); */
  /* printf("LL^T:\n"); */
  /* print_matrix(H); */

  // Clean-up
  gsl_eigen_symmv_free(e);
  gsl_matrix_free(H);

  
  //count up number of non-zero modes.
  // non-zero means having a frequency greater than 10^-8

  modes = 0;

  // the order in the test is important b/c (modes < f->size) will short-circut
  // and prevent GSL out of range errors
  while ((modes < f->size) && (gsl_vector_get(f,modes) > 1e-8)){modes++;}
  debug("Of %d modes, found %d with non-zero frequency.", (int) f->size, (int) modes);


  // Allocate energies, vector gives bounds checking
  gsl_vector * energies; //freed at end
  energies = gsl_vector_calloc(modes);   // energies of the dim-6 (or w/e) modes

  double E_total = 0.0;

  // Assign mode energies
  // modes - 1: last mode assigned at end deterministically
  for (size_t i = 0; i < (modes - 1); i++){
    double Ei = 0.0; //not available outside of block
    double Ri = gsl_rng_uniform_pos(r); //not available outside of block
    // i + 1.0: Hase's formula starts with i = 1 
    Ei = (E0 - E_total)*(1.0 - pow(Ri, 1.0 / ((double) modes - (1.0 + (double)i) ) ) );
    E_total += Ei;
    gsl_vector_set(energies, i, Ei);
  }
  gsl_vector_set(energies, modes - 1, E0 - E_total); // last mode; modes - 1: indexing

  gsl_vector * normal_q, * normal_p;
  normal_q = gsl_vector_calloc(dim);
  normal_p = gsl_vector_calloc(dim);
 
  // Assign Phases to normal coordinates (zero frequency modes are assigned to 0 (really, left at 0))
  for (size_t i = 0; i < modes; i++){
    double Ri = gsl_rng_uniform_pos(r);

    double Ei = gsl_vector_get(energies, i);
    double wi = sqrt(gsl_vector_get(f, i));
    // printf("Ri\t\tEi\twi\n%g,\t%g,\t%g\n", Ri, Ei, wi);
    gsl_vector_set(normal_q, i,
		   (sqrt(2.0*Ei)/wi)*cos(2*PI*Ri));
    gsl_vector_set(normal_p, i,
		   -1.0*sqrt(2.0*Ei)*sin(2*PI*Ri));
  }

  gsl_vector_set_zero(phase_space_point);
  // views for working with the phase-space point
  gsl_vector_view phase_x = gsl_vector_subvector(phase_space_point, (size_t)   0,   dim);
  gsl_vector_view phase_v = gsl_vector_subvector(phase_space_point,          dim,   dim);

  gsl_matrix_set_zero(temp);
  gsl_blas_dgemm((CBLAS_TRANSPOSE_t) CblasNoTrans
		 , (CBLAS_TRANSPOSE_t) CblasNoTrans
		 ,1.0
		 , M2v
		 , L
		 , 0.0
		 , temp);

  gsl_vector_memcpy(&phase_x.vector, x0); //dest, src
  gsl_blas_dgemv((CBLAS_TRANSPOSE_t) CblasNoTrans
		 , 1.0
		 , temp
		 , normal_q
		 , 1.0
		 , &phase_x.vector);

  gsl_blas_dgemv((CBLAS_TRANSPOSE_t) CblasNoTrans
		 , 1.0
		 , temp
		 , normal_p
		 , 0.0
		 , &phase_v.vector);


  E_total=(*v)(&phase_x.vector)+kinetic_energy(&phase_v.vector,mass); //E_total is now the real energy and not the target from the normal modes

  debug("E_total:\t%g",E_total);
  debug("E0:\t\t%g",E0);
  debug("Deviation:\t%g", fabs((E_total-E0)/E0));

 

  assert (dim % 3 == 0); //time to make sure we have a system in 3-space;

  size_t centers = dim / 3; 

  gsl_vector * j_addition, * omega, *tempv, * cm_vel;
  gsl_matrix * inverse_inertia_tensor, * temp3, * temp3a;
  j_addition = gsl_vector_alloc((size_t) 3); //freeed at end
  omega      = gsl_vector_alloc((size_t) 3); //freeed at end
  cm_vel     = gsl_vector_alloc((size_t) 3); //freeed at end
  tempv      = gsl_vector_alloc((size_t) 3); //freeed at end
  inverse_inertia_tensor = gsl_matrix_alloc((size_t) 3, (size_t) 3); //freeed at end
  temp3                  = gsl_matrix_alloc((size_t) 3, (size_t) 3); //freeed at end
  temp3a                 = gsl_matrix_alloc((size_t) 3, (size_t) 3); //freeed at end


  unsigned int counter  = 0;
  do{
    counter++;
    debug("Entering correction round %d.", counter);

    get_angular_momentum(&phase_x.vector, mass, &phase_v.vector, j_addition);
    //printf("  Angular Momentum:\t");
    //print_vect(j_addition);
    gsl_vector_scale(j_addition, -1.0); //since we target 0 angular momentum

    // construct and invert inertia tensor
    get_inertia_tensor(&phase_x.vector, mass, inverse_inertia_tensor);
    invert_3x3(inverse_inertia_tensor);

    // construct angular velocity
    gsl_blas_dgemv((CBLAS_TRANSPOSE_t) CblasNoTrans, 1.0, inverse_inertia_tensor, j_addition, 0.0, omega);

    for (size_t center = 0; center < centers; center++){
      gsl_vector_view vi = gsl_vector_subvector(&phase_v.vector, 3 * center, (size_t) 3);//velocity of each center
      gsl_vector_view ri = gsl_vector_subvector(&phase_x.vector, 3 * center, (size_t) 3);//position of each center
      gsl_vector_memcpy(tempv, omega);
      vector_cross(tempv, &ri.vector); //tempv hold Del-v
      gsl_blas_daxpy(1.0, tempv, &vi.vector);
    }

    E_total=(*v)(&phase_x.vector)+kinetic_energy(&phase_v.vector,mass);
    debug("  Energy Deviation:\t%g", fabs((E_total-E0)/E0));

    //rescale x,v
    for (size_t center = 0; center < centers; center++){
      gsl_vector_view vi   = gsl_vector_subvector(&phase_v.vector, 3 * center, (size_t) 3);//velocity of each center
      gsl_vector_view ri   = gsl_vector_subvector(&phase_x.vector, 3 * center, (size_t) 3);//position of each center
      gsl_vector_view ri_0 = gsl_vector_subvector(x0, 3 * center, (size_t) 3);//position of each center
      gsl_vector_scale(&vi.vector, sqrt(E0/E_total));
      gsl_vector_memcpy(tempv, &ri.vector);
      gsl_vector_sub(tempv, &ri_0.vector);// a <- a-b

      gsl_vector_memcpy(&ri.vector, &ri_0.vector);
      gsl_blas_daxpy(sqrt(E0/E_total),tempv, &ri.vector);
    }
    
    //compute CM motion
    center_of_mass_velocity(&phase_v.vector, mass, cm_vel);
    //printf("  CM velocity:\t\t");
    //print_vect(cm_vel);

    //subtract off CM motion
    for (size_t center = 0; center < centers; center++){
      gsl_vector_view vi = gsl_vector_subvector(&phase_v.vector, 3 * center, (size_t) 3);//velocity of each center
      gsl_vector_sub(&vi.vector, cm_vel);
    }

    E_total=(*v)(&phase_x.vector)+kinetic_energy(&phase_v.vector,mass);
    get_angular_momentum(&phase_x.vector, mass, &phase_v.vector, j_addition);
    center_of_mass_velocity(&phase_v.vector, mass, cm_vel);
    
    if (counter > 1e3){//assume we've hit a fixed point or some other foolishness
      return GSL_FAILURE;
    }

  }while((energy_tolerance           < fabs((E_total-E0)/E0)) ||
	 (angular_momentum_tolerance < gsl_blas_dnrm2(j_addition)) ||
	 (cm_velocity_tolerance      < gsl_blas_dnrm2(cm_vel))
	 );

  //printf("---Correction rounds done---\n");
  debug("Energy Deviation:\t%g", fabs((E_total-E0)/E0));
  //get_angular_momentum(&phase_x.vector, mass, &phase_v.vector, j_addition);
  //printf("Angular Momentum:\t");
  //print_vect(j_addition);
  //printf("CM velocity:\t\t");
  //print_vect(cm_vel);

  
  gsl_vector_free(j_addition);
  gsl_vector_free(omega);
  gsl_vector_free(cm_vel);
  gsl_vector_free(tempv);
  gsl_matrix_free(inverse_inertia_tensor);
  gsl_matrix_free(temp3);
  gsl_matrix_free(temp3a);
  gsl_vector_free(normal_q);
  gsl_vector_free(normal_p);
  gsl_vector_free(energies);
  gsl_vector_free(f);
  gsl_matrix_free(L);
  gsl_matrix_free(temp);
  gsl_matrix_free(M2v);

  return GSL_SUCCESS;
}

void construct_inverse_root_mass_matrix(gsl_vector * mass, gsl_matrix * invrootmass){
  size_t dim;
  dim = mass -> size;
  assert((dim == invrootmass -> size1) && (dim == invrootmass -> size2));

  gsl_vector_view invrootmass_diag = gsl_matrix_diagonal(invrootmass);

  gsl_matrix_set_zero(invrootmass);

  for (size_t i = 0; i < dim; i++){
    gsl_vector_set(&invrootmass_diag.vector, i,
		   1/sqrt(gsl_vector_get(mass, i)));
  }
}


// maybe better to break this up...
/* void get_normal_modes(){ */

/* } */
