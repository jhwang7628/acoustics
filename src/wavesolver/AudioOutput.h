#ifndef AUDIO_OUTPUT_H
#define AUDIO_OUTPUT_H
#include <fstream> 
#include <memory>
#include "TYPES.h" 
#include "config.h"

//##############################################################################
// class AudioOutput
//##############################################################################
class AudioOutput
{
    static std::unique_ptr<AudioOutput> _instance; 
    std::ofstream _stream; 
    FloatArray _buffer; 

    AudioOutput(); 
    
public: 
    static std::unique_ptr<AudioOutput> &instance() {return _instance;}
    void SetBufferSize(const int N) {_buffer.resize(N,0.0);} 
    void ResetBuffer();
    void AccumulateBuffer(const FloatArray &data); 

    // IO
    void OpenStream(const std::string &filename); 
    void CloseStream(); 
    bool WriteAndResetBuffer(); 
}; 
using AudioOutput_UPtr = std::unique_ptr<AudioOutput>; 

#endif
