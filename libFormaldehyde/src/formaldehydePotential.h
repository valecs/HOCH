#ifndef FORMALDEHYDE_POTENTIAL_H
#define FORMALDEHYDE_POTENTIAL_H

#ifndef EXTGETPOT
#define EXTGETPOT getpot_
#endif

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <potentials.h>

extern const potential * const hochPotential;
extern double potential_formaldehyde(const gsl_vector * r);
extern void   gradient_formaldehyde(const gsl_vector * r, gsl_vector * grad);

extern const potential * const hochPotential_nop;
extern double potential_formaldehyde_nop(const gsl_vector * r);
extern void   gradient_formaldehyde_nop(const gsl_vector * r, gsl_vector * grad);

#endif
