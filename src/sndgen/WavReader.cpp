#include <sndgen/WavReader.hpp>

//##############################################################################
// Deleter for SNDFILE used in smarter pointer.
//##############################################################################
void DeleteSNDFILE(SNDFILE *sndFile)
{
    if (sndFile != nullptr)
        sf_close(sndFile); 
}
