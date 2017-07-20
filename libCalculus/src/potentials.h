#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <gsl/gsl_vector.h>

typedef double (*energy)(const gsl_vector *);
typedef void (*gradient)(const gsl_vector *, gsl_vector *);

typedef struct potential {
  energy    v;
  gradient  g;
} potential;

#endif
