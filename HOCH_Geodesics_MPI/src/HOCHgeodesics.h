#ifndef HOCHGEODESICS_H
#define HOCHGEODESICS_H

#include <gsl/gsl_vector.h>
#include <linkedList.h>

llist * getOptimizedHOCHGeodesic(const gsl_vector * A, const gsl_vector * B, double landscapeEnergy);
llist * findHOCHGeodesic(const gsl_vector * A, const gsl_vector * B, double landscapeEnergy);
llist * optimizeHOCHGeodesic(llist * path, double landscapeEnergy);
void path_free(llist * path);
void resetCenterOfMass(gsl_vector * v);

#endif
