#include <modal_model/ImpulseSeriesObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <wavesolver/AccelerationNoiseVibrationalSource.h>
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
}

//##############################################################################
// This function add impulses. The vertex index should be for surface triangle
// mesh (not volumetric tetrahedron).
//##############################################################################
void ImpulseSeriesObject::
AddImpulse(const ImpulseSeriesObject::ImpactRecord &record)
{
    // need the notion of mesh to check whether appliedVertex makes sense
    if (_objectMesh == nullptr)
        throw std::runtime_error("**ERROR** Adding impulse to an object whose mesh is not defined yet.");
    if (record.appliedVertex >= _objectMesh->num_vertices())
        throw std::runtime_error("**ERROR** Adding impulse whose applied vertex ID ("+std::to_string(record.appliedVertex)+") does not make sense (mesh has only "+std::to_string(_objectMesh->num_vertices())+" vertices).");

    _lastImpulseTime = record.timestamp;
    _firstImpulseTime = std::min<REAL>(_firstImpulseTime, record.timestamp);
    _impulses.push_back(record);
}

//##############################################################################
// This function gets single impulse frame by index
//##############################################################################
void ImpulseSeriesObject::
GetImpulse(const int &index, REAL &timestamp, int &vertex, Vector3d &impulse)
{
    const auto &record = _impulses.at(index);
    timestamp = record.timestamp;
    vertex = record.appliedVertex;
    impulse = record.impactVector;
}

//##############################################################################
// This function gets all impulses that have timestamps within [timeStart, timeStop]
//##############################################################################
void ImpulseSeriesObject::
GetImpulse(const REAL &timeStart, const REAL &timeStop, std::vector<ImpactRecord> &records)
{
    if (timeStart > timeStop || timeStart > _lastImpulseTime || timeStop < _firstImpulseTime){
        return;}

    records.clear();
    const int N_impulses = N_Impulses();
    for (int frame_idx=0; frame_idx<N_impulses; ++frame_idx)
    {
        const auto &record = _impulses.at(frame_idx);
        const REAL &timestamp = record.timestamp;
        if (timestamp >= timeStart && timestamp <= timeStop)
            records.push_back(record);
    }
}

//##############################################################################
// This function gets all impulses that have support at this query time point
//##############################################################################
void ImpulseSeriesObject::
GetImpulseWithinSupport(const REAL &timeStart, std::vector<ImpactRecord> &records)
{
    if (timeStart > _lastImpulseTime && timeStart < _firstImpulseTime){
        return;}

    records.clear();
    const int N_impulses = N_Impulses();
    for (int frame_idx=0; frame_idx<N_impulses; ++frame_idx)
    {
        const auto &record = _impulses.at(frame_idx);
        if (timeStart >= record.timestamp && timeStart < record.timestamp+record.supportLength)
            records.push_back(record);
    }
}

//##############################################################################
// This function returns a vector of forces that is nonzero at time 'time' using
// Hertz contact model.
//
// Using Hertz model, we can spread out a given impulse and the force profile is
// assumed to be f(t) = gamma*S(t;t_0,tau)*J_hat, where
//   J_hat is the normalized impulse J
//   gamma = pi norm(J)/(2 tau) is the scaling factor so that
//                              norm(\int_0^{\infty} f(t)dt              )
//                            =      \int_0^{\infty} gamma*S(t;t_0,tau)dt
//                            = norm(J)
//   S(t; t_0, tau) is a half sine that has support in [t_0, t_0+tau]
//##############################################################################
void ImpulseSeriesObject::
GetForces(const REAL &time, std::vector<ImpactRecord> &records)
{
    GetImpulseWithinSupport(time, records);
    const int N = records.size();
    for (int frame_idx=0; frame_idx<N; ++frame_idx)
    {
        auto &r = records.at(frame_idx);
        const REAL j = r.impactVector.length();
        r.impactVector = (r.impactVector/j) *
            r.gamma * AccelerationNoiseVibrationalSource::ComputeS(r,time);
    }
}

//##############################################################################
//##############################################################################
void ImpulseSeriesObject::
GetRangeOfImpulses(REAL &firstImpulseTime, REAL &lastImpulseTime)
{
    firstImpulseTime = _firstImpulseTime;
    lastImpulseTime = _lastImpulseTime;
}

//##############################################################################
// Threshold the impulse magnitude by impact velocity. See 2012 PAN paper.
//##############################################################################
void ImpulseSeriesObject::
Filter()
{
    // in-place erase low magnitude impulses
    for (std::vector<ImpactRecord>::iterator it=_impulses.begin(); it<_impulses.end();)
    {
        if (abs(it->contactSpeed) < IMPULSE_VEL_THRESHOLD)
            it = _impulses.erase(it);
        else
            ++it;
    }

    // update cached field
    if (N_Impulses() > 0)
    {
        _firstImpulseTime = _impulses.at(0).timestamp;
        _lastImpulseTime = (N_Impulses() == 1 ? _firstImpulseTime + SMALL_NUM : _impulses.at(N_Impulses()-1).timestamp);
    }
    else
    {
        _firstImpulseTime = 0.0;
        _lastImpulseTime = 0.0;
    }
}

//##############################################################################
//##############################################################################
REAL ImpulseSeriesObject::
GetFirstImpulseTime(const REAL &startTime) const
{
    if (startTime < _firstImpulseTime || _impulses.size()==0)
        return _firstImpulseTime;
    REAL first = std::numeric_limits<REAL>::max();
    for (const auto &ir : _impulses)
        first = std::min(first, ir.timestamp);
    return first;
}

//##############################################################################
//##############################################################################
void ImpulseSeriesObject::
PrintAllImpulses()
{
    for (const auto &imp : _impulses)
        std::cout << imp << std::endl;
}

//##############################################################################
//##############################################################################
std::ostream &operator <<(std::ostream &os, const ImpulseSeriesObject::ImpactRecord &record)
{
    os << "--------------------------------------------------------------------------------\n"
       << "Struct ImpulseSeriesObject::ImpactRecord\n"
       << "--------------------------------------------------------------------------------\n"
       << " timestamp                           : " << record.timestamp << "\n"
       << " impact vector                       : " << record.impactVector << "\n"
       << " support length (contact time scale) : " << record.supportLength << "\n"
       << " contact speed                       : " << record.contactSpeed << "\n"
       << " applied vertex id                   : " << record.appliedVertex << "\n"
       << "--------------------------------------------------------------------------------"
       << std::flush;
    return os;
}
