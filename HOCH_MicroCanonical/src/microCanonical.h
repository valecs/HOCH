#ifndef MICROCANONICAL_H
#define MICROCANONICAL_H

#include <potentials.h>
#include <gsl/gsl_vector.h>

int get_phase_space_point(gsl_vector  * x0   // initial position
			   ,gsl_vector * mass // vector of masses
			   ,energy       v    // potential
			   ,double       E0   // total desired energy
			   ,gsl_vector * phase_space_point
			   ,double       hessian_tolerance
			   ,double       energy_tolerance
			   ,double       angular_momentum_tolerance
			   ,double       cm_velocity_tolerance);

#endif
