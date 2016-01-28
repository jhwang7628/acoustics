#ifndef MULTIPOLE_MATH_H
#define MULTIPOLE_MATH_H

#include <config.h>
#include <complex>

namespace Multipole{
	std::complex<REAL> spherical_harmonics(int m, int n, REAL theta, REAL phi);
	std::complex<REAL> hankel_1st(int l, REAL z);
	std::complex<REAL> hankel_2nd(int l, REAL z);
	REAL spherical_bessel(int l, REAL z);
	std::complex<REAL> regular_basis(int m, int n, REAL k, REAL r, REAL theta, REAL phi);
	std::complex<REAL> regular_basis_dir_deriv(int m, int n, REAL k, REAL r, REAL theta, REAL phi, REAL dx, REAL dy, REAL dz);
}

#endif