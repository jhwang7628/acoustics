#include <wavesolver/FDTD_RigidObject.h> 
#include <io/TglMeshReader.hpp>
#include <utils/SimpleTimer.h>
#include <utils/Conversions.h>

//##############################################################################
//##############################################################################
void FDTD_RigidObject::
Initialize()
{
    assert(_parsed); 

    // building mesh 
    _mesh.reset(new TriangleMesh<REAL>()); 
    std::cout << _meshScale << std::endl;
    if (MeshObjReader::read(_meshFileName.c_str(), *_mesh, false, false, _meshScale)==SUCC_RETURN)
        _mesh->generate_normals(); 
    else 
        throw std::runtime_error("**ERROR** Cannot read mesh from"+_meshFileName);

    // building signed distance field in the original frame of reference
    _signedDistanceField.reset(
            DistanceFieldBuilder::BuildSignedClosestPointField(
                _meshFileName.c_str(), 
                _signedDistanceFieldResolution, 
                _signedDistanceFieldFilePrefix.c_str()
                )
            );
}

//##############################################################################
// Employ a linear time (in surface element) check to update bounding box
//##############################################################################
void FDTD_RigidObject::
UpdateBoundingBox()
{
    if (_mesh)
    {
        Point3<double> minBound, maxBound; 
        _mesh->bounding_box(minBound, maxBound); 
        _bboxWorld.Update(minBound,maxBound); 
    }
    else 
        throw std::runtime_error("**ERROR** mesh not set");
}

//##############################################################################
//##############################################################################
REAL FDTD_RigidObject::
DistanceToMesh(const double &x, const double &y, const double &z)
{
    if (!_signedDistanceField)
        throw std::runtime_error("**ERROR** distance field not built.");

    if (!_bboxWorld.Inside(x,y,z,1.1))
        return std::numeric_limits<REAL>::lowest();

    Eigen::Vector3d position(x,y,z); 
    position = _modelingTransformInverse*position;
    return _signedDistanceField->distance(position[0],position[1],position[2]); 
}

//##############################################################################
// Note that transforming vector and point has different syntax due to Eigen
// API. Example (vec1 -> vec2) 
//
//  For points: vec2 = transformation          * vec1
//  For vector: vec2 = transformation.linear() * vec1
//
//##############################################################################
bool FDTD_RigidObject::
NormalToMesh(const double &x, const double &y, const double &z, Vector3d &queriedNormal)
{
    if (!_signedDistanceField)
        throw std::runtime_error("**ERROR** distance field not built.");
    if (!_bboxWorld.Inside(x,y,z,1.1))
        return false;
    Eigen::Vector3d position(x,y,z); 
    position = _modelingTransformInverse.linear()*position;
    queriedNormal = _signedDistanceField->gradient(Conversions::ToVector3<double>(position));
    return true; 
}

//##############################################################################
//##############################################################################
REAL FDTD_RigidObject::
EvaluateBoundaryCondition(const Vector3d &boundaryPoint, const Vector3d &boundaryNormal, const REAL &time)
{
    REAL bcValue=0; 
    const SourceIterator sourceEnd = _vibrationalSources.end(); 
    for (SourceIterator it=_vibrationalSources.begin(); it!=sourceEnd; ++it) 
    {
        bcValue += (*it)->Evaluate(boundaryPoint, boundaryNormal, time);
    }
    return bcValue; 
}

//##############################################################################
//##############################################################################
void FDTD_RigidObject::
PrintBoundingBox()
{
    std::cout << "minBound = " << _bboxWorld.minBound << std::endl; 
    std::cout << "maxBound = " << _bboxWorld.maxBound << std::endl; 
}

//##############################################################################
//##############################################################################
void FDTD_RigidObject::
TestQueryDistance()
{
    const int N = 3; 
    const REAL xMin = -0.1; 
    const REAL yMin = -0.1; 
    const REAL cellSize = (-xMin)*2/(double)N;
    for (int ii=0; ii<N; ++ii) 
        for (int jj=0; jj<N; ++jj) 
        {
            const REAL x = xMin + (double)ii*cellSize; 
            const REAL y = yMin + (double)jj*cellSize; 
            std::cout << "distance(" << x << ", " << y << ", 0) = " << DistanceToMesh(x,y,0) << std::endl; 
        }
}

//##############################################################################
//##############################################################################
void FDTD_RigidObject::
TestObjectBoundaryCondition()
{
    const Vector3d boundaryPoint(0,0,0); 
    const Vector3d boundaryNormal(1,1,1); 
    const REAL dt = 1E-6;
    const int N = 20000;
    for (int ii=0; ii<N; ii++)
    {
        std::cout << "---\n";
        COUT_SDUMP((double)ii*dt);
        const REAL result = EvaluateBoundaryCondition(boundaryPoint, boundaryNormal, (double)ii*dt); 
        COUT_SDUMP(result);
    }
}


