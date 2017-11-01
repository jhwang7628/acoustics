#include <map>
#include "wavesolver/FDTD_ShellObject.h"
//##############################################################################
// Static Initialization
//##############################################################################
std::queue<boost::filesystem::path> FDTD_ShellObject::_DataBuffer::qPrefixes; 
const REAL FDTD_ShellObject::_DataStep::stepSize = 1./88200.;

//##############################################################################
// Function DistanceToMesh
//   Always return INF since SDF is not well-defined in ShellObject. If signed 
//   distance to nearest mesh triangles are needed, call 
//   ComputeClosestPointOnMesh in TriangleMeshKDTree classs. I don't want to 
//   relay that call right here because I implicitly assumes DistanceToMesh
//   is an O(1) operation and calls it many times. 
//##############################################################################
REAL FDTD_ShellObject:: 
DistanceToMesh(const Vector3d &position)
{
    if (!_init)
        Initialize();
    return std::numeric_limits<REAL>::max(); 
}

//##############################################################################
// Function ReflectAgainstBoundary
//##############################################################################
int FDTD_ShellObject:: 
ReflectAgainstBoundary(const Vector3d &originalPoint, 
                       Vector3d &reflectedPoint, 
                       Vector3d &boundaryPoint, 
                       Vector3d &erectedNormal, 
                       REAL &distanceTravelled, 
                       const int &startFromTriangle)
{
    throw std::runtime_error("**ERROR** Should not call ReflectAgainstBoundary on Shells"); 
}

//##############################################################################
// Function ReflectAgainstBoundary
//##############################################################################
int FDTD_ShellObject:: 
ReflectAgainstBoundary(const Vector3d &originalPoint, 
                       const std::set<int> &triangles,
                       Vector3d &reflectedPoint, 
                       Vector3d &boundaryPoint, 
                       Vector3d &erectedNormal, 
                       REAL &distance)
{
    int closestTriangle;
    std::vector<REAL> alld; 
    const bool hasDeformation = GetDisplacements(alld); 
    std::vector<int> cands; 
    std::copy(triangles.begin(), triangles.end(), std::back_inserter(cands));
    Vector3d tmp_pp; 
    distance = GetMeshPtr()->ComputeClosestPointOnMeshHelper(originalPoint,
                                                             cands, 
                                                             boundaryPoint,
                                                             closestTriangle, 
                                                             tmp_pp,
                                                             hasDeformation,
                                                             alld); 
    reflectedPoint = boundaryPoint + (originalPoint - boundaryPoint)*2.0; 
    erectedNormal = (originalPoint - boundaryPoint).normalized();
    return closestTriangle; 
}

//##############################################################################
// Function UpdateBoundingBox
//##############################################################################
void FDTD_ShellObject:: 
UpdateBoundingBox()
{
    _bboxWorldUnion2Steps = _bboxWorld; 
    Point3d b_l(std::numeric_limits<REAL>::max(),
                std::numeric_limits<REAL>::max(),
                std::numeric_limits<REAL>::max()); 
    Point3d b_h(std::numeric_limits<REAL>::lowest(),
                std::numeric_limits<REAL>::lowest(),
                std::numeric_limits<REAL>::lowest());
    for (int ii=0; ii<_mesh->num_vertices(); ++ii)
    {
        const Vector3d v = GetVertexPos(ii); 
        b_l.x = std::min(b_l.x, v.x);
        b_l.y = std::min(b_l.y, v.y);
        b_l.z = std::min(b_l.z, v.z);
        b_h.x = std::max(b_h.x, v.x);
        b_h.y = std::max(b_h.y, v.y);
        b_h.z = std::max(b_h.z, v.z);
    }
    for (int ii=0; ii<3; ++ii)
    {
        if (EQUAL_FLOATS(b_l[ii], b_h[ii]))
        {
            b_l[ii] = b_h[ii] - 1E-8;
            b_h[ii] = b_h[ii] + 1E-8;
        }
    }
    Point3d b_c = (b_l + b_h)/2.0;
    const Vector3d h = (b_h - b_l)/2.0*1.05; // enlarge by 105%
    b_l = b_c - h; 
    b_h = b_c + h; 
    _bboxWorld.Update(b_l, b_h); 
    _bboxWorldUnion2Steps.Union(_bboxWorld);
}
       
//##############################################################################
// Function _BuildDistanceField
//##############################################################################
void FDTD_ShellObject::
Initialize()
{
    // ready displacement/acceleration data
    _dataBuffer.ReadMetaData();
    _init = true; 
}
  
//##############################################################################
// Function UpdatePosAcc
//##############################################################################
void FDTD_ShellObject:: 
UpdatePosAcc(const REAL time)
{
    if (!EQUAL_FLOATS(time, _dataBuffer.lastQueriedTime)) 
        _dataBuffer.FindData(time, &_currentData); 
}

//##############################################################################
// Function GetAllVertexPos
//##############################################################################
void FDTD_ShellObject:: 
GetAllVertexPos(std::vector<Vector3d> &allp) const
{
    allp.clear(); 
    allp.resize(_mesh->num_vertices()); 
    for (int ii=0; ii<allp.size(); ++ii) 
        allp.at(ii) = GetVertexPos(ii); 
}

//##############################################################################
// Function GetAllVertexPos
//##############################################################################
bool FDTD_ShellObject:: 
GetDisplacements(std::vector<REAL> &alld) const
{
    alld.resize(_mesh->num_vertices()*3);
    if (_currentData) 
    {
        const auto d = _currentData->displacement; 
        std::copy(d.begin(), d.end(), alld.begin()); 
        return true; 
    }
    std::fill(alld.begin(), alld.end(), (REAL)0.0);
    return false; 
}

//##############################################################################
// Function GetVertexPos
//##############################################################################
Vector3d FDTD_ShellObject:: 
GetVertexPos(const int v_idx) const
{
    Vector3d pos = _mesh->vertex(v_idx); 
    if (_currentData)
    {
        pos.x += _currentData->displacement.at(v_idx*3  ); 
        pos.y += _currentData->displacement.at(v_idx*3+1); 
        pos.z += _currentData->displacement.at(v_idx*3+2); 
    }
    return pos;
}

//##############################################################################
// Function GetVertexAcc
//##############################################################################
Vector3d FDTD_ShellObject:: 
GetVertexAcc(const int v_idx) const 
{
    Vector3d acc; 
    if (_currentData)
    {
        acc.x += _currentData->acceleration.at(v_idx*3  ); 
        acc.y += _currentData->acceleration.at(v_idx*3+1); 
        acc.z += _currentData->acceleration.at(v_idx*3+2); 
    }
    return acc;
}

//##############################################################################
// Function _DataBuffer::FindData
//   Find the data corresponds to time. First scan through all data in current
//   buffer, if cannot find it and time > last frame time in the buffer, read
//   in the next batch. Do so until find the data or no data to read. 
//   
//   @return bool whether data is found
//   @param time required data time
//   @param data output data if found
//##############################################################################
bool FDTD_ShellObject::_DataBuffer::
FindData(const REAL time, _DataStep **data)
{
    assert(time >= bufStartTime); 
    const int offset = (int) std::floor(
            (time - bufStartTime)/FDTD_ShellObject::_DataStep::stepSize); 
    if (offset < bufValidLen) // data in current buffer
    {
        *data = &(buf.at(offset)); 
        lastQueriedTime = time; 
        return true; 
    }
    else if (!qPrefixes.empty()) // read more recursively
    {
        ReadNextBuffer(); 
        return FindData(time, data); 
    }
    else if (bufValidLen > 0) // clamp to last read step data
    {
        *data = &(buf.at(bufValidLen - 1)); 
        lastQueriedTime = time; 
        return true;
    }
    else // no buffer at all, don't return anything
    {
        *data = nullptr; 
        return false; 
    }
}

//##############################################################################
// Function _DataBuffer::ReadNextBuffer
//##############################################################################
void FDTD_ShellObject::_DataBuffer::
ReadNextBuffer()
{
    assert(owner && !owner->_o2iMap.empty() && !owner->_i2oMap.empty());
    bufValidLen = 0;
    for (int bb=0; bb<bufLen && !_DataBuffer::qPrefixes.empty(); ++bb)
    {
        auto p = std::move(_DataBuffer::qPrefixes.front()); 
        const std::string f_dis = p.replace_extension(disSuffix).string(); 
        const std::string f_acc = p.replace_extension(accSuffix).string(); 
        const std::string f_cdis = p.replace_extension(cdisSuffix).string(); 
        const std::string f_cacc = p.replace_extension(caccSuffix).string(); 

        // read displacement and acceleration
        _DataStep &step = buf.at(bb); 
        auto &dis = step.displacement; 
        auto &acc = step.acceleration; 
        auto &cdis = step.constraint_displacement; 
        auto &cacc = step.constraint_acceleration; 
        int N_unconstrained_dof; 
        int N_constrained_dof; 
        bool has_constraint_dis = false; 
        bool has_constraint_acc = false; 

        {
            std::ifstream stream(f_dis.c_str(), std::ios::binary); 
            if (!stream) 
                throw std::runtime_error("**ERROR** cannot open displacement.");
            stream.read((char*)&N_unconstrained_dof, sizeof(int)); 
            dis.resize(N_unconstrained_dof); 
            stream.read((char*)&(dis[0]), sizeof(REAL)*N_unconstrained_dof); 
        }

        {
            std::ifstream stream(f_acc.c_str(), std::ios::binary); 
            if (!stream) 
                throw std::runtime_error("**ERROR** cannot open acceleration.");
            stream.read((char*)&N_unconstrained_dof, sizeof(int)); 
            acc.resize(N_unconstrained_dof); 
            stream.read((char*)&(acc[0]), sizeof(REAL)*N_unconstrained_dof); 
        }

        {
            std::ifstream stream(f_cdis.c_str(), std::ios::binary); 
            if (stream) 
            {
                stream.read((char*)&N_constrained_dof, sizeof(int)); 
                cdis.resize(N_constrained_dof); 
                stream.read((char*)&(cdis[0]), sizeof(REAL)*N_constrained_dof); 
                has_constraint_dis = true; 
            }
        }

        {
            std::ifstream stream(f_cacc.c_str(), std::ios::binary); 
            if (stream) 
            {
                stream.read((char*)&N_constrained_dof, sizeof(int)); 
                cacc.resize(N_constrained_dof); 
                stream.read((char*)&(cacc[0]), sizeof(REAL)*N_constrained_dof); 
                has_constraint_acc = true; 
                assert(has_constraint_dis);
            }
        }

        // reorder vertices based on map
        {
            const int N_total_v         = owner->_i2oMap.size();
            const int N_unconstrained_v = N_unconstrained_dof/3; 
            const int N_constrained_v   = N_total_v - N_unconstrained_v; 
            std::vector<REAL> tmp_d(N_total_v*3, (REAL)0); 
            std::vector<REAL> tmp_a(N_total_v*3, (REAL)0); 
            for (int ii=0; ii<N_unconstrained_v; ++ii)
            {
                const int o_v = owner->_i2oMap.at(ii); 
                tmp_d.at(o_v*3  ) = dis.at(ii*3  );
                tmp_d.at(o_v*3+1) = dis.at(ii*3+1);
                tmp_d.at(o_v*3+2) = dis.at(ii*3+2);
                tmp_a.at(o_v*3  ) = acc.at(ii*3  );
                tmp_a.at(o_v*3+1) = acc.at(ii*3+1);
                tmp_a.at(o_v*3+2) = acc.at(ii*3+2);
            }
            if (has_constraint_dis && has_constraint_acc)
            {
                assert(N_constrained_dof = N_constrained_v*3);
                for (int jj=N_unconstrained_v; jj<N_total_v; ++jj)
                {
                    const int o_v = owner->_i2oMap.at(jj); 
                    const int ii = jj-N_unconstrained_v;
                    tmp_d.at(o_v*3  ) = cdis.at(ii*3  );
                    tmp_d.at(o_v*3+1) = cdis.at(ii*3+1);
                    tmp_d.at(o_v*3+2) = cdis.at(ii*3+2);
                    tmp_a.at(o_v*3  ) = cacc.at(ii*3  );
                    tmp_a.at(o_v*3+1) = cacc.at(ii*3+1);
                    tmp_a.at(o_v*3+2) = cacc.at(ii*3+2);
                }
            }
            std::swap(tmp_d, dis); 
            std::swap(tmp_a, acc); 

            // compute and store the list of constrained vertices
            step.constrainedVertices.resize(N_constrained_v); 
            for (int jj=0; jj<N_constrained_v; ++jj)
            {
                const int ii = N_total_v-1-jj;
                step.constrainedVertices.at(jj) = owner->_i2oMap.at(ii); 
            }
        }
        step.frameID = std::stoi(p.stem().string()); 
        step.frameTime = (REAL)step.frameID * 
                         FDTD_ShellObject::_DataStep::stepSize; 
        ++bufValidLen; 
        _DataBuffer::qPrefixes.pop();
    }
    if (bufValidLen > 0)
        bufStartTime = buf.at(0).frameTime; 
}

//##############################################################################
// Function ReadMetaData
//##############################################################################
void FDTD_ShellObject::_DataBuffer::
ReadMetaData()
{
    using namespace boost::filesystem; 
    path p(owner->_dataDir.c_str()); 
    std::vector<std::string> ordered_filenames;
    std::map<std::string, path> map_path; 
    path p_map;
    for (directory_iterator it(p); it!=directory_iterator(); ++it)
    {
        const auto filename = it->path().filename().string(); 
        if (filename.find("."+disSuffix) != std::string::npos)
        {
            const auto stemstr = it->path().stem().string(); 
            ordered_filenames.push_back(stemstr);
            auto p_tmp = it->path(); 
            map_path[stemstr] = p_tmp.replace_extension();
        }
        else if (filename.find(owner->_vertexMapFile) 
                 != std::string::npos)
        {
            p_map = it->path(); 
        }
    }
    std::sort(ordered_filenames.begin(), ordered_filenames.end(), 
            [](const std::string &a, const std::string &b)
              {return std::stoi(a) < std::stoi(b);});

    // queue all the filenames (prefix with path)
    for (const auto key : ordered_filenames)
        _DataBuffer::qPrefixes.push(map_path.at(key)); 
    std::cout << "Queued shell data step count: " 
              << _DataBuffer::qPrefixes.size() << std::endl;

    // read vertex map 
    {
        std::ifstream stream(p_map.string(), std::ios::in); 
        const int N_vert = owner->_mesh->num_vertices(); 
        if (stream) 
        {
            std::string line; 
            owner->_o2iMap.reserve(N_vert); 
            owner->_i2oMap.reserve(N_vert); 
            int buf[2];
            while (std::getline(stream, line))
            {
                std::istringstream iss(line); 
                iss >> buf[0] >> buf[1]; 
                owner->_o2iMap.push_back(buf[0]); 
                owner->_i2oMap.push_back(buf[1]);
            }
        }
        else // no vertex map found, assume no mapping
        {
            std::cout << "No vertex map found, use identity mapping\n";
            owner->_o2iMap.resize(N_vert); 
            owner->_i2oMap.resize(N_vert); 
            for (int ii=0; ii<N_vert; ++ii)
            {
                owner->_o2iMap.at(ii) = ii; 
                owner->_i2oMap.at(ii) = ii; 
            }
        }
    }
}

