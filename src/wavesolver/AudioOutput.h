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
    std::fstream stream; 
    FloatArray data; 
    AudioOutput() = default; 
    
public: 
    static std::unique_ptr<AudioOutput> &instance()
    {return _instance;}
    // IO
    void OpenStream(const std::string &filename); 
    void CloseStream(); 
}; 
using AudioOutput_UPtr = std::unique_ptr<AudioOutput>; 

#endif
