#include "igl/read_triangle_mesh.h"
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
    // read wav file to set up speaker
    _data.clear();
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

    // read markers for positioning control
    if (_pControl->type == "markers")
    {
        ReadMarkersData();
    }
}

//##############################################################################
// Function ReadMarkersData
//##############################################################################
void SpeakerPressureSource::
ReadMarkersData()
{
    std::cout << "Reading markers for the speaker pressure sources ...\n";
    auto control =
        std::dynamic_pointer_cast<MarkersPositioningControl>(_pControl);
    using namespace boost::filesystem;
    // iterate and find all marker files
    path p(control->dir.c_str());
    std::vector<path> filenames;
    for (directory_iterator it(p); it!=directory_iterator(); ++it)
        if (it->path().extension().string() == ".obj")
            filenames.push_back(it->path());
    auto ParseFileID = [](const path &a)
    {
        std::vector<std::string> tokens;
        boost::split(tokens, a.stem().string(), [](char c){return c == '_';});
        return std::stoi(tokens.at(1));
    };
    std::sort(filenames.begin(), filenames.end(),
            [&](const path &a, const path &b){
                return ParseFileID(a) < ParseFileID(b);});
    control->markersFile = std::move(filenames);
    std::cout << " Done. Total marker files read: "
              << control->markersFile.size()
              << std::endl;
}

//##############################################################################
// Function UpdateBoundingBox
//##############################################################################
void SpeakerPressureSource::
UpdateBoundingBox()
{
    _bboxWorld.Update(_pControl->position - _widthSpace*GAUSSIAN_CHECK_BOUND,
                      _pControl->position + _widthSpace*GAUSSIAN_CHECK_BOUND);
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

    const Vector3d &pos = _pControl->position;
    REAL value = pow((evaluatePosition.x - pos.x),2.0)
               + pow((evaluatePosition.y - pos.y),2.0)
               + pow((evaluatePosition.z - pos.z),2.0);
    value  = -value / (2.0*pow(_widthSpace,2.0));
    value  = exp(value);

    const REAL delayTime = time - (evaluatePosition-pos).length()/_c;
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

//##############################################################################
//##############################################################################
bool SpeakerPressureSource::
UpdateTime(const REAL time)
{
    bool changed = false;
    if (_pControl->type == "markers")
    {
        auto control =
            std::dynamic_pointer_cast<MarkersPositioningControl>(_pControl);
        const int idx =
            (int)((time - control->startTime)*(REAL)control->frameRate);
        if (idx >= 0 && idx < control->markersFile.size())
        {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            const boost::filesystem::path &p = control->markersFile.at(idx);
            igl::read_triangle_mesh(p.string().c_str(), V, F);
            const Eigen::Vector3d C = V.colwise().sum() / (REAL)V.rows();
            control->position.set(C[0], C[1], C[2]);
            changed = true;
        }
    }
    _pControl->time = time;
    return changed;
}
