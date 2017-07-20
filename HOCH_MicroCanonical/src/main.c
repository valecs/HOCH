#include "microCanonical.h"
#include <potentials.h>
#include <formaldehydePotential.h>
#include <formaldehydeProperties.h>
#include "pathHelpers.h"
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <gsl/gsl_vector.h>

#include <math.h>
#include <assert.h>

#define STRATT_ROOT "STRATT_ROOT"
#define DBLF "%#22.16g"
#define BUF 1024

#define HACM 219474.6313705

// Allows setting the starting seed with the energy parameter.
extern unsigned long int gsl_rng_default_seed;

int main (int argc, char ** argv){
  const size_t dim = formaldehydeNumDimensions;
  FILE * f;

  char fname[BUF];

  const char * stratt_root = getenv(STRATT_ROOT);

  if (!stratt_root){
    fprintf(stderr, "Variable, " STRATT_ROOT ", not found; exiting!\n");
    exit(-1);
  }

  int E0, total;

  if (argc != 3){
    fprintf(stderr, "Requires 2 arguments: ENERGY COUNT\n");
    fprintf(stderr, "ENERGY is in units of wavenumbers.\n");
    exit(3);
  }

  //in cm^-1 , converted to Ha in call to get_phase_space_point
  E0 = atoi(argv[1]);
  total = atoi(argv[2]);
  
  fprintf(stderr, "Generating %d points at an energy of %d wavenumbers.\n", total, E0);

  gsl_rng_default_seed = (long unsigned int) E0;
  
  if (snprintf(fname, BUF, "%s/data/HOCH/properties/"
	       "hoch-EQ.geo.vector", stratt_root) >= BUF ){
    fprintf(stderr, "Output truncated; exiting!\n");
    exit (1);
  }

  gsl_vector * x0;
  x0 = gsl_vector_alloc(dim);
  f = fopen(fname, "r");
  if (!f){
    fprintf(stderr, "Unable to open equilibrium geometry, %s; exiting!\n", fname);
    exit(1);
  }
  else{
    gsl_vector_fscanf(f,x0);
    fclose (f);
  }

  if (snprintf(fname, BUF, "%s/data/HOCH/properties/"
	       "hoch.mass.vector", stratt_root) >= BUF ){
    fprintf(stderr, "Output truncated; exiting!\n");
    exit (1);
  }
  
  gsl_vector * mass;
  mass = gsl_vector_alloc(dim);
  f = fopen(fname, "r");
  if (!f){
    fprintf(stderr, "Unable to open mass vector, %s; exiting!\n", fname);
    exit(1);
  }
  else{
    gsl_vector_fscanf(f,mass);
    fclose (f);
  }

  energy v = hochPotential->v;
  
  gsl_vector * phase_space_point;
  phase_space_point = gsl_vector_alloc(dim*2);

  double hessian_tolerance          = 1e-4; // in units of a0
  double energy_tolerance           = 1e-8; // in fractional deviation
  double angular_momentum_tolerance = 1e-8; // max mag of L vector
  double cm_velocity_tolerance      = 1e-8; // max mag of R vector

  for (int counter = 0; counter < total; counter++){
    if (GSL_SUCCESS == get_phase_space_point(x0, mass, v, ((double) E0/HACM)
					     , phase_space_point
					     , hessian_tolerance
					     , energy_tolerance
					     , angular_momentum_tolerance
					     , cm_velocity_tolerance)){

      if (snprintf(fname, BUF, "%s/data/HOCH/geometries/"
		   "%05d-%05d.gamma", stratt_root, E0, counter + 1) >= BUF ){
	fprintf(stderr, "Output truncated; exiting!\n");
	exit (1);
      }

      printf("Phase space point found; writing to: %s.\n", fname);
    
      f = fopen(fname, "w");
      if (!f){
	fprintf(stderr, "Unable to write to vector, %s; exiting!\n", fname);
	exit(2);
      }
      else{
	gsl_vector_fprintf(f, phase_space_point, DBLF);
	fclose(f);
      }
    }
    else{
      printf("No phase-space point found for round %d; trying again!\n", counter + 1);
      counter--;
    }
  }

  return 0;
}

