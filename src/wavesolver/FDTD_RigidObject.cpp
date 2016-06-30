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
    UpdateBoundingBox();
}

//##############################################################################
// Employ a linear time (in surface element) check to update bounding box
// union of the last two updates are also computed and stored
//##############################################################################
void FDTD_RigidObject::
UpdateBoundingBox()
{
    // cache the old one
    _bboxWorldUnion2Steps = _bboxWorld; 
    // compute the new one
    if (_mesh)
    {
        Vector3<REAL> minBound(D_INF, D_INF, D_INF); 
        Vector3<REAL> maxBound(-D_INF, -D_INF, -D_INF); 
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
        _bboxWorld.Update(minBound, maxBound); 
    }
    else 
        throw std::runtime_error("**ERROR** mesh not set");
    // union step
    _bboxWorldUnion2Steps.Union(_bboxWorld);
}

//##############################################################################
//##############################################################################
REAL FDTD_RigidObject::
DistanceToMesh(const double &x, const double &y, const double &z)
{
    if (!_signedDistanceField)
        throw std::runtime_error("**ERROR** distance field not built.");
    if (!_bboxWorld.Inside(x, y, z, AABB_CHECK_TOLERANCE_SCALE))
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
EvaluateBoundaryAcceleration(const Vector3d &boundaryPoint, const Vector3d &boundaryNormal, const REAL &time)
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
REAL FDTD_RigidObject::
EvaluateBoundaryVelocity(const Vector3d &boundaryPoint, const Vector3d &boundaryNormal, const REAL &time)
{
    REAL bcValue=0; 
    const SourceIterator sourceEnd = _vibrationalSources.end(); 
    for (SourceIterator it=_vibrationalSources.begin(); it!=sourceEnd; ++it) 
    {
        bcValue += (*it)->EvaluateVelocity(boundaryPoint, boundaryNormal, time);
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
    assert(_signedDistanceField!=nullptr && (DistanceToMesh(originalPoint.x,originalPoint.y,originalPoint.z)<DISTANCE_TOLERANCE));

    //erectedNormal = _signedDistanceField->gradient(originalPoint); 
    NormalToMesh(originalPoint.x, originalPoint.y, originalPoint.z, erectedNormal);
    erectedNormal.normalize(); 
    distanceTravelled = 2.0 * fabs(DistanceToMesh(originalPoint.x, originalPoint.y, originalPoint.z)); 
    boundaryPoint  = originalPoint + erectedNormal * (0.5*distanceTravelled);
    reflectedPoint = originalPoint + erectedNormal * (    distanceTravelled); 
    const REAL newDistance = DistanceToMesh(reflectedPoint);
    const bool reflectSuccess = (newDistance >= DISTANCE_TOLERANCE); 

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
            std::cerr << "**ERROR** reflected point " << originalPoint << "->" << reflectedPoint << " still inside object : " << newDistance << std::endl; 
    }

    return reflectSuccess;
}

//##############################################################################
// Without a better name, this function finds the image point (IP) for the fresh 
// cell in ghost-cell implementation by extending the current query point (CU)
// to find IP. Boundary point (BI) will also be located for boundary condition
// evaluation later. The fresh cell is contracted to be outside the boundary 
// when this function is called.
//
//           /
// IP  CU  BI
//  o---o---o  inside boundary
//          |  
//          |
//
//##############################################################################
bool FDTD_RigidObject::
FindImageFreshCell(const Vector3d &currentPoint, Vector3d &imagePoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &distanceTravelled)
{
    assert(_signedDistanceField != nullptr);
    //assert(_signedDistanceField!=nullptr && (DistanceToMesh(currentPoint.x,currentPoint.y,currentPoint.z)>DISTANCE_TOLERANCE));

    //erectedNormal = _signedDistanceField->gradient(currentPoint); 
    NormalToMesh(currentPoint.x, currentPoint.y, currentPoint.z, erectedNormal);
    erectedNormal.normalize(); 
    //distanceTravelled = -_signedDistanceField->distance(currentPoint)*2; 
    distanceTravelled = DistanceToMesh(currentPoint.x, currentPoint.y, currentPoint.z);
    if (distanceTravelled > DISTANCE_TOLERANCE) // located correctly outside the boundary
    {
        boundaryPoint = currentPoint - erectedNormal * (distanceTravelled);
        imagePoint    = currentPoint + erectedNormal * (distanceTravelled); 
    }
    else // inside the buffer, follow a similar procedure as the ghost cell method 
    {
        boundaryPoint = currentPoint - erectedNormal * (distanceTravelled);
        imagePoint    = currentPoint - erectedNormal * (2.0*distanceTravelled); // want it to be outside the boundary
    }

    const REAL newDistance = DistanceToMesh(imagePoint);
    const bool isExterior = (newDistance > DISTANCE_TOLERANCE); 
    if (true) 
    {
        Vector3d boundaryNormal;
        NormalToMesh(boundaryPoint.x, boundaryPoint.y, boundaryPoint.z, boundaryNormal); 
        //Vector3d boundaryNormal = _signedDistanceField->gradient(boundaryPoint); 
        boundaryNormal.normalize(); 
        if (erectedNormal.dotProduct(boundaryNormal) < 0.5)
        {
            std::cerr << "**WARNING** erected normal and true normal deviates. This might cause inaccuracy for the imposed Neumann boundary condition at original point : " 
                      << currentPoint 
                      << "; the dot product is : " << erectedNormal.dotProduct(boundaryNormal) << std::endl; 
        }
        if (!isExterior)
            std::cerr << "**ERROR** reflected point " << currentPoint << "->" << imagePoint << " still inside object : " << newDistance << std::endl; 
    }

    return isExterior;
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
        const REAL result = EvaluateBoundaryAcceleration(boundaryPoint, boundaryNormal, (double)ii*dt); 
        std::cout << "result at time " << (double)ii*dt << " is " << result << std::endl;
    }
}

