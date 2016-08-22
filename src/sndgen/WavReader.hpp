#ifndef SOUND_READER_HPP
#define SOUND_READER_HPP

#include <vector>
#include <string>
#include <memory>
#include <sndfile.h>
#include <assert.h>
#include <utils/STL_Wrapper.h>

//##############################################################################
// This class reads data from a single channel WAV file. 
//##############################################################################
template<typename T>
class WavReader
{
    private: 
        SF_INFO                     _sfInfo; 
        std::shared_ptr<SNDFILE>    _sndFile; 

        template<typename TRead> 
            struct type {};
        template<typename TRead>
            sf_count_t ReadAux(std::vector<T> &data, const int &frames);
        sf_count_t ReadAux(type<double>, std::vector<T> &data, const int &frames);
        sf_count_t ReadAux(type<float>, std::vector<T> &data, const int &frames);
        sf_count_t ReadAux(type<short>, std::vector<T> &data, const int &frames);
        sf_count_t ReadAux(type<int>, std::vector<T> &data, const int &frames);

    public: 
        bool Open(const std::string &wavFile);
        void Close();
        void Read(std::vector<T> &data); 
};

//##############################################################################
// Class helper function declaration
//##############################################################################
void DeleteSNDFILE(SNDFILE *sndFile); 

//##############################################################################
//##############################################################################
template<typename T>
bool WavReader<T>::
Open(const std::string &wavFile)
{
    // API: set format to zero before sf_open(), except for the case of RAW file.
    _sfInfo.format = 0;
    _sndFile.reset(sf_open(wavFile.c_str(), SFM_READ, &_sfInfo), DeleteSNDFILE);
    if (!_sndFile)
        throw std::runtime_error("**ERROR** Cannot open wavfile:" + wavFile); 
    
    // print information
    std::cout << "Wav file:            " << wavFile << "\n"
              << " Sampling rate:      " << _sfInfo.samplerate << "\n"
              << " Channels:           " << _sfInfo.channels << "\n"
              << " Frames:             " << _sfInfo.frames << "\n";
    std::cout << std::flush;

    return (_sndFile != nullptr); 
}

//##############################################################################
//##############################################################################
template<typename T>
void WavReader<T>::
Close()
{
    if (_sndFile != nullptr)
        sf_close(_sndFile.get());
}

//##############################################################################
//##############################################################################
template<typename T>
void WavReader<T>::
Read(std::vector<T> &data)
{
    assert(_sndFile); 
    data.resize(_sfInfo.frames * _sfInfo.channels); 

    const sf_count_t frameCount = ReadAux<T>(data, _sfInfo.frames); 
    std::cout << " Actual frames read: " << frameCount << std::endl;
}

//#############################################################################
// This function is a specialized helper function that interfaces sndfile using 
// overload. Ref:
//   http://stackoverflow.com/questions/5512910/explicit-specialization-of-template-class-member-function
//#############################################################################
template<typename T> template<typename TRead>
sf_count_t WavReader<T>::
ReadAux(std::vector<T> &data, const int &frames)
{
    return ReadAux(type<TRead>(), data, frames);
}

//##############################################################################
//##############################################################################
template<typename T>
sf_count_t WavReader<T>::
ReadAux(type<double>, std::vector<T> &data, const int &frames)
{
    return sf_read_double(_sndFile.get(), &data[0], frames) ;
}

//##############################################################################
//##############################################################################
template<typename T>
sf_count_t WavReader<T>::
ReadAux(type<float>, std::vector<T> &data, const int &frames)
{
    return sf_read_float(_sndFile.get(), &data[0], frames) ;
}

//##############################################################################
//##############################################################################
template<typename T>
sf_count_t WavReader<T>::
ReadAux(type<short>, std::vector<T> &data, const int &frames)
{
    return sf_read_short(_sndFile.get(), &data[0], frames) ;
}

//##############################################################################
//##############################################################################
template<typename T>
sf_count_t WavReader<T>::
ReadAux(type<int>, std::vector<T> &data, const int &frames)
{
    return sf_read_int(_sndFile.get(), &data[0], frames) ;
}

#endif
