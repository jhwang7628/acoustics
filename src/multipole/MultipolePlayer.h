#ifndef MULTIPOLE_PLAYER_H
#define MULTIPOLE_PLAYER_H

#include <sndgen/player/SoundPlayer.h>
#include "MultipoleData.h"

#include <cmath>
#include <complex>
#include <vector>

	class MultipolePlayer : public SoundPlayer{
	public:
		MultipolePlayer(const MultipoleData & data,
			double sampleRate=44100, double framesPerBuffer=10000)
			: SoundPlayer(sampleRate, framesPerBuffer), multipoleData(data){

		}
	
		virtual float soundAtTime(float t);
		void listenAt(double x, double y, double z);
	public:
		MultipoleData multipoleData;
		std::map<int, std::pair<std::complex<double>, double> > estimates;
	};

#endif