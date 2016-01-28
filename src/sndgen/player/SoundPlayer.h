#ifndef SOUNDPLAYER_H 
#define SOUNDPLAYER_H

	class SoundPlayerImpl;
	
	class SoundPlayer{
	public:
		static void wait(double seconds);
		SoundPlayer(double sampleRate=44100, double framesPerBuffer=256);
		~SoundPlayer();

		void start(float t=0);
		void stop();
		virtual float soundAtTime(float t);
	private:
		double _sampleRate;
		double _framesPerBuffer;
		double _time;

		SoundPlayerImpl * impl;
	};

#endif