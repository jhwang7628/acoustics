#include "wavesolver/SpeakerPressureSource.h"
#include "wavesolver/SpeakerVibrationalSource.h"
#include "sndgen/WavReader.hpp"

//#define USE_DDT_WAV_AS_SRC

//##############################################################################
// Function Initialize
//##############################################################################
void SpeakerPressureSource::
Initialize(const std::string &speakerFile)
{
    _data.clear();
    // read wav file to set up speaker
    WavReader<REAL> reader;
    try
    {
        reader.Open(speakerFile);
    }
    catch (std::runtime_error &e)
    {
        std::cout << "**WARNING** Cannot open file for SpeakerPressureSource: "
                  << speakerFile << std::endl;
        return;
    }
#ifdef USE_DDT_WAV_AS_SRC
    std::vector<REAL> p;
    reader.ReadChannel(p, 0);  // only read one channel
    // compute derivative
    _data.resize(p.size()-1);
    for (int ii=0; ii<_data.size(); ++ii)
        _data.at(ii) = (p.at(ii+1) - p.at(ii))/(1./_sampleRate);
#else
    reader.ReadChannel(_data, 0);
#endif
    _sampleRate = reader.SampleRate();
    STL_Wrapper::PrintVectorContent(std::cout, _data, 10);
    reader.Close();
    _dataFile = speakerFile;
}

//##############################################################################
// Function UpdateBoundingBox
//##############################################################################
void SpeakerPressureSource::
UpdateBoundingBox()
{
    _bboxWorld.Update(_position - _widthSpace*GAUSSIAN_CHECK_BOUND,
                      _position + _widthSpace*GAUSSIAN_CHECK_BOUND);
}

//##############################################################################
// Use bounding box to check if evaulation is needed
//##############################################################################
REAL SpeakerPressureSource::
Evaluate(const Vector3d &evaluatePosition, const Vector3d &normal,
         const REAL &time, const int &hintTriangle)
{
    if (!_bboxWorld.Inside(evaluatePosition))
        return 0.0;

    REAL value = pow((evaluatePosition.x - _position.x),2.0)
               + pow((evaluatePosition.y - _position.y),2.0)
               + pow((evaluatePosition.z - _position.z),2.0);
    value  = -value / (2.0*pow(_widthSpace,2.0));
    value  = exp(value);

    const REAL delayTime = time - (evaluatePosition-_position).length()/_c;
    const int idx = (delayTime - _startTime)/_sampleRate;
    if (idx >= 0 && idx < _data.size())
        value *= _data.at(idx);
    else
        value = 0.0;

    return value;
}

//##############################################################################
//##############################################################################
REAL SpeakerPressureSource::
Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time)
{
    throw std::runtime_error("**ERROR** Cannot sample SpeakerPressureSource using vertexID");
    return 0.0;
}

//##############################################################################
//##############################################################################
Vector3d SpeakerPressureSource::
Evaluate(const int &vertexID, const REAL &time)
{
    throw std::runtime_error("**ERROR** Cannot sample SpeakerPressureSource using vertexID");
    return Vector3d();
}

//##############################################################################
//##############################################################################
REAL SpeakerPressureSource::
EarliestEventTime(const REAL startTime) const //TODO
{
    return 0.0;
}

//##############################################################################
//##############################################################################
bool SpeakerPressureSource::
IsZero(const REAL t,
       const bool checkBound,
       const Vector3d &minBound,
       const Vector3d &maxBound) const
{

    if (t < _startTime ||
        t > _data.size()*_sampleRate)
        return true;
    const int idx = (int)((t - _startTime)/_sampleRate);
    return fabs(_data.at(idx)) < _dataThreshold;
    return true;
}
