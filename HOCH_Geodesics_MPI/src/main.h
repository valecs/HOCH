#ifndef HOCH_Geodesics_MPI_H
#define HOCH_Geodesics_MPI_H

#include <gsl/gsl_vector.h>
#include "kinds.h"

struct pointsHolder{
  gsl_vector * A;
  gsl_vector * B;
  gsl_vector * C;
};

struct params{
  int manifestID;
  geodesicKind kind;
  int E0;
  const char * stratt_root;
};


#endif
