#include "wavesolver/AudioOutput.h" 
//##############################################################################
// Static member definition
//##############################################################################
AudioOutput_UPtr AudioOutput::_instance = AudioOutput_UPtr(new AudioOutput); 

//##############################################################################
// Constructor
//##############################################################################
AudioOutput::AudioOutput()
{
}

//##############################################################################
// Function ResetBuffer
//##############################################################################
void AudioOutput::
ResetBuffer()
{
    std::fill(_buffer.begin(), _buffer.end(), 0.0); 
}

//##############################################################################
// Function AccumulateBuffer
//##############################################################################
void AudioOutput::
AccumulateBuffer(const FloatArray &data)
{
    assert(data.size() == _buffer.size()); 
    std::transform(data.begin(), data.end(), _buffer.begin(), _buffer.begin(),
                   std::plus<REAL>()); 
}

//##############################################################################
// Function OpenStream
//##############################################################################
void AudioOutput::
OpenStream(const std::string &filename)
{
    _stream.open(filename.c_str(), std::ios::binary);
}

//##############################################################################
// Function CloseStream
//##############################################################################
void AudioOutput::
CloseStream()
{
    _stream.close();
}

//##############################################################################
// Function WriteAndResetBuffer
//##############################################################################
bool AudioOutput::
WriteAndResetBuffer()
{
    if (!_stream.is_open())
        return false; 
    _stream.write((char*) &(_buffer[0]), sizeof(REAL)*_buffer.size()); 
    _stream.flush();

    COUT_SDUMP(_buffer[0]); 
    ResetBuffer(); 
    return true;
}
