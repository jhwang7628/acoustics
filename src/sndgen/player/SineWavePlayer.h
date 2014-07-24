#ifndef SINEWAVEPLAYER_H
#define SINEWAVEPLAYER_H

#include "SoundPlayer.h"

#include <cmath>

	class SineWavePlayer : public SoundPlayer{
	public:
		SineWavePlayer(double freq, double amplitude, double phase,
					   double sampleRate=44100, double framesPerBuffer=256)
			: SoundPlayer(sampleRate, framesPerBuffer){
			this->_freq = 2*acos(-1)*freq;
			this->_amp = amplitude;
			this->_phase = phase;
		}
	
		virtual float soundAtTime(float t){
			return 0.0001*_amp*sin(_freq*t + _phase);
		}

		void setFrequency(double freq){
			this->_freq = 2*acos(-1)*freq;
		}

		void setAmplitude(double amp){
			this->_amp = amp;
		}

		void setPhase(double phase){
			this->_phase = phase;
		}
	private:
		double _freq;
		double _amp;
		double _phase;
	};

#endif