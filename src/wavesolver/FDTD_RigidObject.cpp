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
        Point3<REAL> minBound(D_INF, D_INF, D_INF); 
        Point3<REAL> maxBound(-D_INF, -D_INF, -D_INF); 
        Eigen::Vector3d pointBuffer; 
        std::vector<Point3<REAL>> &meshVertices = _mesh->vertices(); 
        const typename std::vector<Point3<REAL>>::const_iterator end = meshVertices.end();
        for(typename std::vector<Point3<REAL>>::const_iterator it=meshVertices.begin(); it!=end; ++it)
        {
            pointBuffer[0] = it->x; 
            pointBuffer[1] = it->y; 
            pointBuffer[2] = it->z; 
            pointBuffer = _modelingTransform * pointBuffer; 

            minBound.x = min(minBound.x, pointBuffer[0]);
            minBound.y = min(minBound.y, pointBuffer[1]);
            minBound.z = min(minBound.z, pointBuffer[2]);

            maxBound.x = max(maxBound.x, pointBuffer[0]);
            maxBound.y = max(maxBound.y, pointBuffer[1]);
            maxBound.z = max(maxBound.z, pointBuffer[2]);
        }
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
        return std::numeric_limits<REAL>::max();

    Eigen::Vector3d position(x,y,z); 
    position = _modelingTransformInverse*position;
    return _signedDistanceField->distance(position[0],position[1],position[2]); 
}

//##############################################################################
//##############################################################################
REAL FDTD_RigidObject::
DistanceToMesh(const Vector3d &position)
{
    return DistanceToMesh(position.x, position.y, position.z); 
}

//##############################################################################
// Note that transforming vector and point has different syntax due to Eigen
// API. Example (vec1 -> vec2) 
//
//  For points: vec2 = transformation          * vec1
//  For vector: vec2 = transformation.linear() * vec1
//
//  TODO 
//  too many conversions in this call
//
//##############################################################################
bool FDTD_RigidObject::
NormalToMesh(const double &x, const double &y, const double &z, Vector3d &queriedNormal)
{
    if (!_signedDistanceField)
        throw std::runtime_error("**ERROR** distance field not built.");
    if (!_bboxWorld.Inside(x,y,z,1.1))
    {
        const REAL limit = std::numeric_limits<REAL>::max(); 
        queriedNormal = Vector3d(limit, limit, limit); 
        return false;
    }

    Eigen::Vector3d position(x,y,z); 
    position = _modelingTransformInverse*position;
    queriedNormal = _signedDistanceField->gradient(Conversions::ToVector3<double>(position));
    Eigen::Vector3d normal = Conversions::ToEigen<double>(queriedNormal); 
    normal = _modelingTransform.linear()*normal; 
    queriedNormal = Conversions::ToVector3(normal); 
    return true; 
}

//##############################################################################
//##############################################################################
bool FDTD_RigidObject::
NormalToMesh(const Vector3d &position, Vector3d &queriedNormal)
{
    return NormalToMesh(position.x, position.y, position.z, queriedNormal);
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
// Reflect the given point against the boundary. check if the reflected point 
// is indeed outside the boundary (of current object). 
//
// produce warning if the following scenario occurs 
//  1. erected normal deviates from the boundary normal by too much. 
//  2. image point is still inside the geometry. 
//##############################################################################
bool FDTD_RigidObject::
ReflectAgainstBoundary(const Vector3d &originalPoint, Vector3d &reflectedPoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &distanceTravelled)
{
    assert(_signedDistanceField!=nullptr && (DistanceToMesh(originalPoint.x,originalPoint.y,originalPoint.z)<0));

    //erectedNormal = _signedDistanceField->gradient(originalPoint); 
    NormalToMesh(originalPoint.x, originalPoint.y, originalPoint.z, erectedNormal);
    erectedNormal.normalize(); 
    //distanceTravelled = -_signedDistanceField->distance(originalPoint)*2; 
    distanceTravelled = -1.0 * DistanceToMesh(originalPoint.x, originalPoint.y, originalPoint.z); 
    boundaryPoint  = originalPoint + erectedNormal * (0.5*distanceTravelled);
    reflectedPoint = originalPoint + erectedNormal * (    distanceTravelled); 
    const bool reflectSuccess = (_signedDistanceField->distance(reflectedPoint) > DISTANCE_TOLERANCE); 

    if (true) 
    {
        Vector3d boundaryNormal;
        NormalToMesh(boundaryPoint.x, boundaryPoint.y, boundaryPoint.z, boundaryNormal); 
        //Vector3d boundaryNormal = _signedDistanceField->gradient(boundaryPoint); 
        boundaryNormal.normalize(); 
        if (erectedNormal.dotProduct(boundaryNormal) < 0.5)
        {
            std::cerr << "**WARNING** erected normal and true normal deviates. This might cause inaccuracy for the imposed Neumann boundary condition at original point : " 
                      << originalPoint 
                      << "; the dot product is : " << erectedNormal.dotProduct(boundaryNormal) << std::endl; 
        }
        if (!reflectSuccess)
            std::cerr << "**ERROR** reflected point still inside object" << std::endl; 
    }

    return reflectSuccess;
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
        const REAL result = EvaluateBoundaryCondition(boundaryPoint, boundaryNormal, (double)ii*dt); 
        std::cout << "result at time " << (double)ii*dt << " is " << result << std::endl;
    }
}


