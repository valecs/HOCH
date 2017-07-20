#include <potentials.h>
#include <formaldehydePotential.h>
#include <gradient.h>
#include "debug_print.h"
#include "cmHelpers.h"
#include "hcohHelpers.h"
#include "pathHelpers.h"

#include "mpigrid.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>

#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include <unistd.h>
#include <getopt.h>


#define STRATT_ROOT "STRATT_ROOT"
#define HACM 219474.6313705
#define BUF 1024

int odefunc_hcoh(double t, const double x[], double f[], void * params);
double dist_H2_CO(gsl_vector * x, gsl_vector * mass);
void usage(char * cmd);

int main(int argc, char ** argv){
  const energy potential_energy = hochPotential -> v;
  static const size_t dim = 12;
  FILE * f; 

  const char * stratt_root = getenv(STRATT_ROOT);

  if (!stratt_root){
    fprintf(stderr, "Variable, " STRATT_ROOT ", not found; exiting!\n");
    exit(1);
  }

  if (5 != argc){
    usage(argv[0]);
    exit(1);
  }
  
  int trajectory = getTaskID(argc, argv);
  
  int E0 = -1;
  
  int c = 1;
  opterr = 0;  //since we do our own error-handling
  const char optstring[]="+e:h";
  while (-1 != (c = getopt(argc, argv, optstring))){
    switch (c){
    case 'e':
      E0 = atoi(optarg);
      break;
    case 'h':
      usage(argv[0]);
      exit(0);
      break;
    case '?':
      fprintf(stderr, "Bad argument; exiting!\n");
      exit(1);
      break;
    default:
      fprintf(stderr, "Internal error; exiting!");
      exit(1);
      break;
    }
  }

  char fname[BUF];

  gsl_vector * mass = gsl_vector_calloc (dim);
  // the diagonal of inverse_mass_matrix is 1/mass
  gsl_matrix * inverse_mass_matrix = gsl_matrix_calloc (dim, dim);
    
  {// initialize mass representations
    if (snprintf(fname, (size_t) BUF, "%s/data/HOCH/properties/"
		 "hoch.mass.vector", stratt_root) >= BUF ){
      fprintf(stderr, "Output truncated; exiting!\n");
      exit(2);
    }

    f = fopen(fname, "r");
    if (!f){
      fprintf(stderr, "Unable to open mass vector, %s; exiting!\n", fname);
      exit(3);
    }
    else{
      gsl_vector_fscanf(f, mass);
      fclose (f);
    }
    
    gsl_matrix_set_identity (inverse_mass_matrix);
    gsl_vector_view midiag = gsl_matrix_diagonal(inverse_mass_matrix);
    gsl_vector_div(&midiag.vector, mass);
  }

  gsl_odeiv2_system sys = {
    .function    = &odefunc_hcoh
    , .jacobian  = (int (*) (double , const double*, double * dfdy, double*, void *)) NULL
    , .dimension = 2 * dim
    , .params    = (void *) inverse_mass_matrix
  };

  //double * y = calloc(2 * dim, sizeof(double));
  double y[2 * dim];
  gsl_vector_view y_view = gsl_vector_view_array(y      , dim*2);
  gsl_vector_view x      = gsl_vector_view_array(y      , dim);
  gsl_vector_view v      = gsl_vector_view_array(y + dim, dim);

  const gsl_odeiv2_step_type * stepperT = gsl_odeiv2_step_rk4;

  //these 3 objects are freed at the end of main
  gsl_odeiv2_evolve  * e    = gsl_odeiv2_evolve_alloc(dim * 2);
  gsl_odeiv2_control * con  = gsl_odeiv2_control_y_new(1e-6, 0.0); //epsilon_absolute, epsilon_relative
  gsl_odeiv2_step    * step = gsl_odeiv2_step_alloc (stepperT,dim*2);





  /* //DEBUGGING VARS */
  /* double db_E, db_lmag; */
  /* gsl_vector * db_L = gsl_vector_alloc((size_t) 3); */

  if (snprintf(fname, (size_t) BUF, "%s/data/HOCH/geometries/"
	       "%05d-%05d.gamma", stratt_root, E0, trajectory) >= BUF ){
    fprintf(stderr, "Output truncated; exiting!\n");
    exit(2);
  }

  f = fopen(fname, "r");
  if (!f){
    fprintf(stderr, "Unable to read to vector, %s; exiting!\n", fname);
    exit(3);
  }
  else{
    gsl_vector_fscanf(f, &y_view.vector);
    fclose(f);
  }
  
  /* sprintf(fname,"initial/%05d-%05d.gamma", E0, trajectory); */
  /* debug("loading %s", fname); */
  /* f = fopen(fname, "r"); assert(f != NULL); */
  /* gsl_vector_fscanf(f,&y_view.vector); */
  /* fclose (f); */

  if (snprintf(fname, (size_t) BUF, "%s/data/HOCH/trajectories/%05d/"
	       "%05d.timeseries", stratt_root, E0, trajectory) >= BUF ){
    fprintf(stderr, "Output truncated; exiting!\n");
    exit(2);
  }

  f = fopen(fname, "w");
  if (!f){
    fprintf(stderr, "Unable to write to timeseries, %s; exiting!\n", fname);
    exit(3);
  }
  else{
    fclose(f);
  }
  
  double time = 0.0;
  double h = 5.0;
  bool state = true;
  do{
    int status = gsl_odeiv2_evolve_apply_fixed_step(e, con, step, &sys, &time, h, y);

    if (GSL_SUCCESS == status){
      if (fmod(time, 100.0)  < h ){// write out state;
	debug("TIME %g: writing to %s", time, fname);
	f = fopen(fname, "a");
	if (!f){
	  fprintf(stderr, "Unable to write to timeseries, %s; exiting!\n", fname);
	  exit(3);
	}
	else{
	  fprintf(f,"%g\n",time);
	  gsl_vector_fprintf(f, &y_view.vector, "%22.16g");
	  fclose(f);
	}
      
	//check to see if we should stop
	if (fmod(time, 1000.0) < h ){ //want to go to 12.1ps in .121 fs steps => 1e5 steps
	  if (time > 1e5)                        {state = false; printf("%05d: TIME\n",trajectory); continue;}
	  if (dist_H2_CO(&x.vector, mass) > 12.0){state = false; printf("%05d: DIST\n",trajectory); continue;}
	}
      }
      if (h < 5.0){
	h /= 0.8;
      }
    }
    else{
      if (h < 1.0){
	printf ("%05d: error, h=%g\n", trajectory, h);
	break;
      }
      else{
	h *= 0.8;
	continue;
      }
    }

  }while(state);

  (void) v;
  (void) potential_energy;
  (void) y_view;

  gsl_odeiv2_step_free(step);
  gsl_odeiv2_control_free(con);
  gsl_odeiv2_evolve_free(e);

  //free(y);

  gsl_matrix_free(inverse_mass_matrix);
  gsl_vector_free(mass);
  
  return 0;
}

int odefunc_hcoh(double t, const double x[], double f[], void * params){
  static const size_t dim = 12;
  energy potfunc_hcoh = hochPotential -> v;

  gsl_matrix * m_inv = (gsl_matrix *)params;

  gsl_vector_const_view r = gsl_vector_const_view_array(x, dim);
  gsl_vector_const_view v = gsl_vector_const_view_array(x + dim, dim);
  gsl_vector_view rdot = gsl_vector_view_array(f, dim);
  gsl_vector_view vdot = gsl_vector_view_array(f + dim, dim);

  gsl_vector_memcpy(&rdot.vector, &v.vector);         //rdot <- v
  grad_5pt(potfunc_hcoh, &r.vector, 1e-4, &vdot.vector);  //vdot <- grad U|r
  gsl_blas_dtrmv((CBLAS_UPLO_t)       CblasUpper
		 ,(CBLAS_TRANSPOSE_t) CblasNoTrans
		 ,(CBLAS_DIAG_t)      CblasNonUnit
		 , m_inv, &vdot.vector); //vdot <- m_inv vdot
  gsl_vector_scale(&vdot.vector, -1.0); //vdot <- -vdot

  (void) t;

  return GSL_SUCCESS;
}

void usage(char * cmd){
  printf("Usage: %s -e ENERGY -t n-m \n",cmd);
  printf("  ENERGY :   energy of trajectory in wavenumbers.\n");
  printf("  n-m:       range of initial conditions to generate trajectories over; [n,m].\n");
  printf("\n");
  printf("Initial conditions for ENERGY and RUN must exist in ./initial \n");
}
