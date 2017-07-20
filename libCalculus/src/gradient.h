#ifndef GRADIENTS_H
#define GRADIENTS_H

#include "potentials.h"
#include <gsl/gsl_vector.h>

int grad_5pt(energy f, const gsl_vector * x0, double h, gsl_vector * grad);

#endif
