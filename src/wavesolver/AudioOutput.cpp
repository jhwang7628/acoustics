#include "wavesolver/AudioOutput.h" 
//##############################################################################
// Static member definition
//##############################################################################
AudioOutput_UPtr AudioOutput::_instance = AudioOutput_UPtr(new AudioOutput); 

//##############################################################################
// Function OpenStream
//##############################################################################
void AudioOutput::
OpenStream(const std::string &filename)
{
    stream.open(filename.c_str(), std::ios::binary|std::ios::out);
}

//##############################################################################
// Function CloseStream
//##############################################################################
void AudioOutput::
CloseStream()
{
    stream.close();
}
