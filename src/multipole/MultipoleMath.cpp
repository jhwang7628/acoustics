#include "MultipoleMath.h"

#include <gsl/gsl_sf.h>

namespace Multipole{
	std::complex<REAL> spherical_harmonics(int m, int n, REAL theta, REAL phi){
		return (m & 1) ? REAL(-gsl_sf_legendre_sphPlm(n, abs(m), cos(theta))) * std::complex<REAL>(cos(m*phi), sin(m*phi))
					   : REAL( gsl_sf_legendre_sphPlm(n, abs(m), cos(theta))) * std::complex<REAL>(cos(m*phi), sin(m*phi));
	}

	REAL spherical_bessel(int l, REAL z){
		return REAL(gsl_sf_bessel_jl(l, z));
	}

	REAL spherical_neumann(int l, REAL z){
		return REAL(gsl_sf_bessel_yl(l, z));
	}

	std::complex<REAL> hankel_1st(int l, REAL z){
		return std::complex<REAL>(spherical_bessel(l, z), spherical_neumann(l, z));
	}

	std::complex<REAL> hankel_2nd(int l, REAL z){
		return std::complex<REAL>(spherical_bessel(l, z), -spherical_neumann(l, z));
	}

	std::complex<REAL> regular_basis(int m, int n, REAL k, REAL r, REAL theta, REAL phi){
		if ( n < abs(m) ) return std::complex<REAL>(0, 0);
		return spherical_bessel(n, k*r) * spherical_harmonics(m, n, theta, phi);
	}

	REAL Amnd(int m, int n){
		int mm = abs(m);
		return n < mm ? 0 : (sqrt(REAL((mm+n+1)*(n-mm+1)) / REAL((2*n+1)*(2*n+3))));
	}

	REAL Bmnd(int m, int n){
		return n < abs(m) ? 0 :
			(m >= 0 ? (sqrt(REAL((n-m-1)*(n-m)) / REAL((2*n-1)*(2*n+1))))
					:-(sqrt(REAL((n-m-1)*(n-m)) / REAL((2*n-1)*(2*n+1))))
					);
	}

	std::complex<REAL> regular_basis_dir_deriv(int m, int n, REAL k, 
		REAL r, REAL theta, REAL phi, REAL dx, REAL dy, REAL dz){
		using namespace std;
		const complex<REAL> A =
			Bmnd(-m-1, n+1) * regular_basis(m+1, n+1, k, r, theta, phi) - 
			Bmnd(m,    n)   * regular_basis(m+1, n-1, k, r, theta, phi);
		const complex<REAL> B =
			Bmnd(m-1, n+1)  * regular_basis(m-1, n+1, k, r, theta, phi) -
			Bmnd(-m,  n)    * regular_basis(m-1, n-1, k, r, theta, phi);

		complex<REAL> rx = (0.5*A + 0.5*B);
		complex<REAL> ry = complex<REAL>((A.imag() - B.imag()) * 0.5, (B.real() - A.real()) * 0.5);
		complex<REAL> rz =
			Amnd(m, n-1) * regular_basis(m, n-1, k, r, theta, phi) -
			Amnd(m, n)   * regular_basis(m, n+1, k, r, theta, phi);
		return (rx*REAL(dx) + ry*REAL(dy) + rz*REAL(dz))*REAL(k);
	}
}