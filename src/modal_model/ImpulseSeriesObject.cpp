#include <modal_model/ImpulseSeriesObject.h> 
//##############################################################################
//
//##############################################################################
ImpulseSeriesObject::
ImpulseSeriesObject()
    : _objectID(-1), 
      _objectMesh(nullptr),
      _firstImpulseTime(std::numeric_limits<REAL>::infinity()), 
      _lastImpulseTime(-1.0)
{
    Initialize();
}

//##############################################################################
//##############################################################################
ImpulseSeriesObject::
ImpulseSeriesObject(const int &objectID, const TriangleMeshPtr &meshPtr)
    : _objectID(objectID), 
      _objectMesh(meshPtr),
      _firstImpulseTime(std::numeric_limits<REAL>::infinity()), 
      _lastImpulseTime(-1.0)
{
    Initialize();
}

//##############################################################################
//##############################################################################
void ImpulseSeriesObject::
Initialize()
{
    _lengthImpulses = 0;
}

//##############################################################################
//##############################################################################
void ImpulseSeriesObject::
AddImpulse(const REAL &timestamp, const int &appliedVertex, const Vector3d &impulse)
{
    _impulseTimestamps.push_back(timestamp); 
    _impulseAppliedVertex.push_back(appliedVertex);
    _impulses.push_back(impulse); 
    if (_impulseTimestamps.size() != _impulseAppliedVertex.size() || 
        _impulseTimestamps.size() != _impulses.size()) 
        throw std::runtime_error("**ERROR** Number of impulse timestamps not equal to number of applied verticesor impulses."); 
    if (timestamp < _lastImpulseTime)  // fine if equal 
        throw std::runtime_error("**ERROR** Added impulse was older than the latest impulse in the series."); 
    else 
        _lastImpulseTime = timestamp;
    _firstImpulseTime = std::min<REAL>(_firstImpulseTime, timestamp); 
    _lengthImpulses = _impulseTimestamps.size(); 
}

//##############################################################################
//##############################################################################
void ImpulseSeriesObject::
GetImpulse(const int &index, REAL &timestamp, int &vertex, Vector3d &impulse)
{
    timestamp = _impulseTimestamps.at(index); 
    vertex = _impulseAppliedVertex.at(index); 
    impulse = _impulses.at(index); 
}

//##############################################################################
//##############################################################################
void ImpulseSeriesObject::
GetImpulse(const REAL &time) // TODO!
{
}

//##############################################################################
//##############################################################################
void ImpulseSeriesObject::
GetImpulseRange(REAL &firstImpulseTime, REAL &lastImpulseTime)
{
    firstImpulseTime = _firstImpulseTime; 
    lastImpulseTime = _lastImpulseTime; 
}
