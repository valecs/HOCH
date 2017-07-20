#ifndef HCOHHELPERS_H
#define HCOHHELPERS_H

#include <gsl/gsl_vector.h>

double dist_H2_CO(gsl_vector * x, gsl_vector * mass);

void get_product_vectors(gsl_vector * g, gsl_vector * m
			 ,gsl_vector * g_HH, gsl_vector * m_HH
			 ,gsl_vector * g_CO, gsl_vector * m_CO);

void get_hoch_center_distances(gsl_vector * x, gsl_vector * r);

#endif
