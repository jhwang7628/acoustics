#include "MultipoleData.h"

#include <cmath>
#include "MultipoleMath.h"
#include <cstdio>

std::complex<REAL> MultipoleData::estimateModeAt(int mode, REAL x, REAL y, REAL z){
	using Multipole::hankel_2nd;
	using Multipole::spherical_harmonics;

	REAL dx = x-cx, dy = y-cy, dz = z-cz;

	REAL r = sqrt(dx*dx + dy*dy + dz*dz);
	REAL theta = acos(dz/r);
	REAL phi = atan2(dy, dx);

	if(this->modes.find(mode) == this->modes.end()){
		return std::complex<REAL>(0, 0);
	}

	MultipoleModeData & data = this->modes[mode];
	int nL = data.numCoefficients();
	
	REAL k = 2*acos(-1.0)*data.frequency()/this->c;

	std::complex<REAL> estimate(0, 0);
	for(int l = 0; l <= nL; l++){
		std::complex<REAL> hank = hankel_2nd(l, r*k);
		for(int m = -l; m <= l; m++){
			estimate += data.coefficient(m, l)*hank*spherical_harmonics(m, l, theta, phi);
		}
	}
	return estimate;
}