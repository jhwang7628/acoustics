#include <boost/algorithm/string.hpp>
#include <sndgen/WavReader.hpp>
#include <wavesolver/SpeakerVibrationalSource.h>
//##############################################################################
//##############################################################################
SpeakerVibrationalSource::DataStep SpeakerVibrationalSource::
ReadObjSeqMetaData(const std::string &dir, const std::string &objPrefix,
                   const std::string &speakerVIdsDir,
                   const std::string &speakerVIdsSuf)
{
    std::cout << "Reading meta data for speaker vibrational source \n";

    using namespace boost::filesystem;

    // helper to parse the filename
    auto ParseFileID = [](const path &a) {
        std::vector<std::string> tokens;
        boost::split(tokens, a.filename().string(), [](char c){return c == '_';});
        return std::stoi(tokens.at(2));};

    // iterate through directory and find main character obj series
    path p(dir.c_str());
    std::vector<path> filenames;
    for (directory_iterator it(p); it!=directory_iterator(); ++it)
    {
        const auto filename = it->path().filename().string();
        if (filename.find(objPrefix) != std::string::npos &&
            it->path().extension().string() == ".obj")
            filenames.push_back(it->path());
    }
    std::sort(filenames.begin(), filenames.end(),
            [&](const path &a, const path &b) {
                const int ia = ParseFileID(a);
                const int ib = ParseFileID(b);
                return ia < ib;
            });

    // find corresponding selected vertices and push data into data steps
    for (auto f : filenames)
    {
        IntArray handles;
        {
            const std::string handlesFile = speakerVIdsDir + "/" + f.stem().string() + speakerVIdsSuf;
            std::ifstream stream(handlesFile.c_str());
            if (stream)
            {
                std::string line;
                int buf;
                while(std::getline(stream, line))
                {
                    std::istringstream iss(line);
                    while (iss >> buf)
                        handles.push_back(buf);
                }
            }
            else
            {
                std::cerr << "**WARNING** Cannot open handles file: "+handlesFile << std::endl;
            }
        }

        DataStep data;
        data.frame = ParseFileID(f);
        data.objFilePrefix = f.stem().string();
        data.handles = std::move(handles);
        _objSeqData.push(std::move(data));
    }

    assert(!_objSeqData.empty());
    DataStep firstStep = std::move(_objSeqData.front());
    assert(firstStep.frame == 0);
    _objSeqData.pop();
    return firstStep;
}

//##############################################################################
//##############################################################################
void SpeakerVibrationalSource::
Initialize(const std::string &speakerFile, const std::vector<int> &handleVIds)
{
    _speakerData.clear();
    _handles.clear();

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
    _speakerWavFile = speakerFile;

    // create set of handles
    _handles.insert(handleVIds.begin(), handleVIds.end());
    PrintHandlesTotalArea(); // FIXME debug
}

//##############################################################################
//##############################################################################
bool SpeakerVibrationalSource::
UpdateTime(const REAL time)
{
    if (_objSeqData.empty())
        return false;

    bool changed = false;
    DataStep step = _objSeqData.front();
    while (!_objSeqData.empty() && time >=
            (REAL)_objSeqData.front().frame*_objSeqSampleRate)
    {
        step = _objSeqData.front();
        _objSeqData.pop();
        changed = true;
    }

    if (changed)
    {
        std::cout << "SpeakerVibrationalSource::Mesh changed: "
                  << step.objFilePrefix << std::endl;
        _owner->Reinitialize(step.objFilePrefix, false, false);
        Initialize(_speakerWavFile, step.handles);
    }

    return changed;
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
    const Vector3d a = Evaluate(vertexID, time);
    return a.dotProduct(normal);
}

//##############################################################################
//##############################################################################
Vector3d SpeakerVibrationalSource::
Evaluate(const int &vertexID, const REAL &time)
{
    REAL a = 0.0;
    if (_handles.find(vertexID) != _handles.end())
    {
        const int idx = (time - _speakerDataStartTime)/_speakerDataSampleRate;
        if (idx >= 0 && idx < _speakerData.size())
            a = _speakerData.at(idx);
    }
    const Vector3d n = _owner->ObjectToWorldVector(
            _owner->GetMeshPtr()->normal(vertexID).normalized());
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

//##############################################################################
//##############################################################################
REAL SpeakerVibrationalSource::
EarliestEventTime(const REAL startTime) const
{
    int i = (startTime - _speakerDataStartTime)/_speakerDataSampleRate;
    while(i<_speakerData.size() &&
          fabs(_speakerData.at(i))<_speakerDataThreshold)
        ++i;
    return _speakerDataStartTime + (REAL)i * _speakerDataSampleRate;
}
//##############################################################################
//##############################################################################
bool SpeakerVibrationalSource::
IsZero(const REAL t,
       const bool checkBound,
       const Vector3d &minBound,
       const Vector3d &maxBound) const
{
    if (t < _speakerDataStartTime ||
        t > _speakerData.size()*_speakerDataSampleRate)
        return true;
    const int idx = (int)((t - _speakerDataStartTime)/_speakerDataSampleRate);
    return fabs(_speakerData.at(idx)) < _speakerDataThreshold;
}

//##############################################################################
//##############################################################################
void SpeakerVibrationalSource::
PrintHandlesTotalArea()
{
    if (!_owner) return;
    auto mesh = _owner->GetMeshPtr();
    auto area = mesh->vertex_areas();
    REAL sum = 0.0;
    for (auto h : _handles)
        sum += area[h];
    std::cout << "SpeakerVibrationalSource::PrintHandlesTotalArea()\n";
    std::cout << " Num handle vertices: " << _handles.size() << "\n";
    std::cout << " Total handle area:   " << sum << "\n";
}
