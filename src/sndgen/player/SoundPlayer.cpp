#include "SoundPlayer.h"
#include <portaudio.h>
#include <cmath>
#include <cstdio>

static int paImplCallback( const void *inputBuffer, void *outputBuffer,
						   unsigned long framesPerBuffer,
						   const PaStreamCallbackTimeInfo* timeInfo,
						   PaStreamCallbackFlags statusFlags,
						   void *userData );

class SoundPlayerImpl
{
public:
	static int numInstances;

	static void checkError(PaError err){
		printf(  "PortAudio error: %s\n", Pa_GetErrorText( err ) );
	}

	static void wait(double seconds){
		Pa_Sleep(seconds*1000);
	}

	SoundPlayerImpl(double sampleRate, double framesPerBuffer, SoundPlayer * parent){
		this->_sampleRate = sampleRate;
		this->_framesPerBuffer = _framesPerBuffer;
		this->_time = 0;
		this->_dt = 1/sampleRate;

		SoundPlayerImpl::numInstances++;
		if(numInstances == 1){
			checkError(Pa_Initialize());
		}
		checkError(Pa_OpenDefaultStream(&stream,
										0,
										1,
										paFloat32,
										sampleRate,
										framesPerBuffer,
										paImplCallback,
										this));
		playing = false;
		_parent = parent;
	}

	~SoundPlayerImpl(){
		if(playing) stop();
		checkError(Pa_CloseStream(stream));

		SoundPlayerImpl::numInstances--;
		if(numInstances == 0){
			checkError(Pa_Terminate());
		}
	}

	void start(double t){
		this->_time = t;
		if(!playing) checkError(Pa_StartStream(this->stream));
		playing = true;
		printf("PLAY\n");
	}

	void stop(){
		if(playing) checkError(Pa_StopStream(this->stream));
		playing = false;
	}

	SoundPlayer * _parent;
	bool playing;
	PaStream * stream;
	double _sampleRate;
	double _framesPerBuffer;
	double _time;
	double _dt;
};

int SoundPlayerImpl::numInstances = 0;

SoundPlayer::SoundPlayer(double sampleRate, double framesPerBuffer){
	this->impl = new SoundPlayerImpl(sampleRate, framesPerBuffer, this);
}

SoundPlayer::~SoundPlayer(){
	delete this->impl;
}

void SoundPlayer::wait(double seconds){
	SoundPlayerImpl::wait(seconds);
}

void SoundPlayer::start(float t){
	this->impl->start(t);
}
void SoundPlayer::stop(){
	this->impl->stop();
}

//Simple sinodal wave
float SoundPlayer::soundAtTime(float t){
	return sin(2*acos(-1)*440*t);
}

static int paImplCallback( const void *inputBuffer, void *outputBuffer,
						   unsigned long framesPerBuffer,
						   const PaStreamCallbackTimeInfo* timeInfo,
						   PaStreamCallbackFlags statusFlags,
						   void *userData ){
	SoundPlayerImpl * s = ((SoundPlayerImpl *) userData);
	SoundPlayer * sp = s->_parent;
	float dt = s->_dt;

	(void) inputBuffer;

	float *out = (float*)outputBuffer;
	for(int i = 0; i < framesPerBuffer; i++){
		*out++ = sp->soundAtTime(s->_time);
		s->_time += dt;
	}
	return 0;
}