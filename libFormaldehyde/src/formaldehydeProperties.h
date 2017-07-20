#ifndef FORMALDEHYDE_PROPERTIES_H
#define FORMALDEHYDE_PROPERTIES_H

#include <stddef.h>

struct offset_t{
  const size_t H1;  const size_t O;  const size_t C;  const size_t H2;
};

struct mass_t{
  const double H1;  const double O;  const double C;  const double H2;
};

extern const struct offset_t formaldehydeOffset;
extern const struct mass_t formaldehydeMass;

extern const double formaldehydeMassArray[4];
extern const size_t formaldehydeOffsetArray[4];

extern const size_t formaldehydeNumDimensions;
extern const int formaldehydeNumCenters;
extern const double formaldehydeMassTotal;

/* use like:
gsl_matrix_const_view M = gsl_matrix_const_view_array(formaldehydeMassMatrixArray, formaldehydeNumDimensions, formaldehydeNumDimensions);
*/
extern const double formaldehydeMassMatrixArray[144];
extern const double formaldehydeInvMassMatrixArray[144];
extern const double formaldehydeCMInvariantProjectionArray[144];

#endif
