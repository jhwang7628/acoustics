#include "MultipolePlayer.h"

#include <utility>
#include <iostream>

void MultipolePlayer::listenAt(double x, double y, double z){
	for(std::map<int,MultipoleData::MultipoleModeData>::iterator it = multipoleData.modes.begin(); it != multipoleData.modes.end(); it++){
		int mode = it->first;
		// std::cout << mode << std::endl;
		double frequency = it->second.frequency();
		// std::cout << frequency << std::endl;
		std::pair<std::complex<double>, double> estimated(multipoleData.estimateModeAt(mode, x, y, z), frequency);
		std::cout << "f: " << frequency << "a: " << abs(estimated.first) << std::endl;
		// std::cout << estimated.first << std::endl;
		this->estimates[mode] = estimated;
	}
	std::cout << std::endl;
}

float MultipolePlayer::soundAtTime(float t){
	double sound = 0;
	double PI = acos(-1);
	double m = 2700;
	double alpha = 6;
	double beta = 1e-7;
	double q = 1;
	for(std::map<int, std::pair<std::complex<double>, double> >::iterator it = estimates.begin(); it != estimates.end(); it++){
		std::pair<std::complex<double>, double> pair = it->second;
		double amplitude = std::abs(pair.first);
		double phase = std::arg(pair.first);
		double frequency = pair.second;
		double omega = 2*PI*frequency;
		double epsi = ((alpha/omega) + beta*omega)/2;
		double epp = sqrt(1-epsi*epsi);

		//sound += amplitude*sin(2*PI*frequency*t + phase);
		sound += amplitude * exp(-epsi*omega*t) *sin(omega*epp*t) * q/(m*omega*epp);
	}
	//sound *= 0.0001*exp(-t*5);
	//std::cout << sound << std::endl;
	return sound;
}