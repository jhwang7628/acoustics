#include <sndgen/WavReader.hpp>
#include <wavesolver/SpeakerVibrationalSource.h> 

//##############################################################################
//##############################################################################
void SpeakerVibrationalSource::
Initialize(const std::string &speakerFile, const std::vector<int> &handleVIds)
{
    // read wav file to set up speaker
    WavReader<REAL> reader; 
    try
    {
        reader.Open(speakerFile); 
    } 
    catch (std::runtime_error &e) 
    {
        std::cout << "**WARNING** Cannot open file for SpeakerVibrationalSource: " 
                  << speakerFile << std::endl; 
        return; 
    }
    reader.ReadChannel(_speakerData, 0);  // only read one channel
    STL_Wrapper::PrintVectorContent(std::cout, _speakerData, 10);
    _speakerDataSampleRate = reader.SampleRate(); 
    reader.Close(); 

    // create set of handles
    _handles.insert(handleVIds.begin(), handleVIds.end());

    // TODO
}

//##############################################################################
//##############################################################################
REAL SpeakerVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle)
{
    // Should not be called
    throw std::runtime_error("**ERROR** this is not implemented");
    return 0.0;
}

//##############################################################################
//##############################################################################
REAL SpeakerVibrationalSource::
Evaluate(const int &vertexID, const Vector3d &normal, const REAL &time)
{
    if (_handles.find(vertexID) != _handles.end())
    {
        const int idx = time/_speakerDataSampleRate; 
        if (idx < 0 || idx >= _speakerData.size())
            return 0.0; 
        else 
            return _speakerData.at(idx);
    }
    return 0.0;
}

//##############################################################################
//##############################################################################
Vector3d SpeakerVibrationalSource::
Evaluate(const int &vertexID, const REAL &time)
{
    const REAL a = Evaluate(vertexID, Vector3d(), time); 
    const Vector3d n = _owner->GetMeshPtr()->normal(vertexID).normalized(); 
    return n*a;
}

//##############################################################################
//##############################################################################
REAL SpeakerVibrationalSource::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented"); 
    return 0.0;
}

//##############################################################################
//##############################################################################
REAL SpeakerVibrationalSource::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented"); 
    return 0.0;
}
