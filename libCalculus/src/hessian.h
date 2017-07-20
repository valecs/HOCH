#ifndef HESSIAN_H
#define HESSIAN_H



#include <gsl/gsl_matrix.h>

void compute_hessian(double (*v)(const gsl_vector *), gsl_matrix * H, const gsl_vector * x0, double delta);
void hessian4       (double (*F)(const gsl_vector *), gsl_matrix * H, const gsl_vector * x0, const double h);

#endif
