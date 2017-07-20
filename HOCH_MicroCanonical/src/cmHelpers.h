#ifndef CMHELPERS_H
#define CMHELPERS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double kinetic_energy(gsl_vector * v, gsl_vector * m);
void get_inertia_tensor(gsl_vector * r, gsl_vector * m, gsl_matrix * I);
void get_angular_momentum(const gsl_vector * r, const gsl_vector * m, const gsl_vector * v, gsl_vector * L);
void invert_3x3(gsl_matrix * M);
void center_of_mass(const gsl_vector * r, const gsl_vector * m, gsl_vector * o);
void center_of_mass_velocity(const gsl_vector * v, const gsl_vector * m, gsl_vector * o);
void vector_cross(gsl_vector * a, gsl_vector * b);
double ext_gsl_vector_dist(gsl_vector * x1, gsl_vector * x2);

#endif
