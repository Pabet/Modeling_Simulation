/*
 * MaxwellBoltzmannDistribution.cpp
 *
 * @Date: 01.10.2010
 * @Author: eckhardw
 */

#include "MaxwellBoltzmannDistribution.h"
#include <cstdlib>

/**
 * helper function for MaxwellBoltzmannDistribution().
 * Generates a gauss deviate, i.e. values according to the normal distribution.
 */
static double GaussDeviate();

void MaxwellBoltzmannDistribution(Particle& p, double factor, int dimensions) {
	auto& v = p.getV();
	for (int i = 0; i < dimensions; i++) {
		v[i] = v[i] + factor * GaussDeviate();
	}
}

// code taken from
//Griebel et. al.: Numerical Simulation in Molecular Dynamics, p. 427
static double GaussDeviate() {
	double a1, a2, s, r, b1;
	static bool iset = false;
	static double b2;

	if (!iset) {
		do {
			a1 = 2.0 * rand() / (RAND_MAX + 1.0) - 1.0;
			a2 = 2.0 * rand() / (RAND_MAX + 1.0) - 1.0;
			r = a1 * a1 + a2 * a2;
		} while (r >= 1.0);
		s = sqrt(-2.0 * log(r) / r);
		b1 = a1 * s;
		b2 = a2 * s;
		iset = true;
		return b1;
	} else {
		iset = false;
		return b2;
	}
}

