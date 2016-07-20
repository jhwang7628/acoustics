#include <modal_model/ImpulseSeriesObject.h> 
//##############################################################################
//
//##############################################################################
ImpulseSeriesObject::
ImpulseSeriesObject()
    : _objectMesh(nullptr),
      _firstImpulseTime(std::numeric_limits<REAL>::infinity()), 
      _lastImpulseTime(-1.0)
{
    Initialize();
}

//##############################################################################
//##############################################################################
ImpulseSeriesObject::
ImpulseSeriesObject(const TriangleMeshPtr &meshPtr)
    : _objectMesh(meshPtr),
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
    // need the notion of mesh to check whether appliedVertex makes sense
    if (_objectMesh == nullptr)
        throw std::runtime_error("**ERROR** Adding impulse to an object whose mesh is not defined yet."); 
    if (appliedVertex >= _objectMesh->num_vertices())
        throw std::runtime_error("**ERROR** Adding impulse whose applied vertex ID ("+std::to_string(appliedVertex)+") does not make sense (mesh has only "+std::to_string(_objectMesh->num_vertices())+" vertices)."); 

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
GetImpulse(const REAL &timeStart, const REAL &timeStop, std::vector<ImpactRecord> &records)
{
    if (timeStart > timeStop || timeStart > _lastImpulseTime || timeStop < _firstImpulseTime){
        return;}

    records.clear(); 
    for (int frame_idx=0; frame_idx<Size(); ++frame_idx) 
    {
        const REAL &timestamp = _impulseTimestamps.at(frame_idx); 
        if (timestamp >= timeStart && timestamp <= timeStop)
        {
            ImpactRecord record; 
            record.impactVector = _impulses.at(frame_idx); 
            record.timestamp = _impulseTimestamps.at(frame_idx); 
            record.appliedVertex = _impulseAppliedVertex.at(frame_idx); 
            records.push_back(record); 
        }
    }
}

//##############################################################################
//##############################################################################
void ImpulseSeriesObject::
GetForces(const REAL &timeStart, const REAL &timeStop, std::vector<ImpactRecord> &records)
{
    GetImpulse(timeStart, timeStop, records); 
    const int N = records.size();
    for (int frame_idx=0; frame_idx<N; ++frame_idx) 
        records.at(frame_idx).impactVector = ConvertImpulseToForce(records.at(frame_idx).impactVector); 
}

//##############################################################################
//##############################################################################
void ImpulseSeriesObject::
GetImpulseRange(REAL &firstImpulseTime, REAL &lastImpulseTime)
{
    firstImpulseTime = _firstImpulseTime; 
    lastImpulseTime = _lastImpulseTime; 
}
