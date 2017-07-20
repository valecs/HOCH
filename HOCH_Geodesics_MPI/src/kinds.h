#ifndef KINDS_H
#define KINDS_H

#include <gsl/gsl_vector.h>

typedef enum {
  nothing,
  direct,
  roaming,
  radical
} geodesicKind;

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
