#include <wavesolver/FDTD_RigidObject.h> 
#include <wavesolver/AccelerationNoiseVibrationalSource.h>
#include <io/TglMeshReader.hpp>
#include <io/TglMeshWriter.hpp>
#include <io/TetMeshReader.hpp>
#include <utils/SimpleTimer.h>
#include <utils/Conversions.h>


//##############################################################################
//##############################################################################
void FDTD_RigidObject::
Initialize(const bool &buildFromTetMesh)
{
    assert(_parsed); 

    if (buildFromTetMesh)
    {
        // first established paths for tet mesh, surface mesh, and correponding sdf
        const std::string tetMeshFile = IO::AssembleFilePath(_workingDirectory, _objectPrefix) + ".tet"; 
        const std::string geoFile = IO::AssembleFilePath(_workingDirectory, _objectPrefix) + ".geo.txt"; 
        const std::string tetSurfaceMeshFile = IO::AssembleFilePath(_workingDirectory, _objectPrefix) + ".tet.obj"; 
        _signedDistanceFieldFilePrefix = tetSurfaceMeshFile + "." + std::to_string(_signedDistanceFieldResolution) + ".dist";

        if (!IO::ExistFile(tetMeshFile))
            throw std::runtime_error("**ERROR** Tet mesh file not exist: " + tetMeshFile); 

        // build/load mesh 
        std::shared_ptr<FixVtxTetMesh<REAL> > tetMesh = std::make_shared<FixVtxTetMesh<REAL> >(); 
        if (FV_TetMeshLoader_Double::load_mesh(tetMeshFile.c_str(), *tetMesh) == SUCC_RETURN) 
        {
            _mesh = std::make_shared<TriangleMeshKDTree<REAL> >(); 
            //_mesh.reset(new TriangleMesh<REAL>()); 
            tetMesh->extract_surface(_mesh.get()); 
            _mesh->generate_normals(); 
            _mesh->generate_mean_curvatures(); 
            _mesh->update_vertex_areas(); 
            _tetMeshIndexToSurfaceMesh = std::make_shared<TetMeshIndexToSurfaceMesh>(); 
            _tetMeshIndexToSurfaceMesh->ReadFromGeoFile(geoFile); 
            if (_tetMeshIndexToSurfaceMesh->N_surfaceVertices() != _mesh->num_vertices())
                throw std::runtime_error("**ERROR** geo file has different number of surface vertices than the surface mesh from tet mesh");
            else 
                std::cout << " Surface mesh and Tet-Surface mapping built.\n";
            _volume = tetMesh->total_volume(); // cache volume
        }
        else 
        {
            throw std::runtime_error("**ERROR** Cannot read mesh from" + tetMeshFile);
        }

        // write the surface mesh
        if (!IO::ExistFile(tetSurfaceMeshFile))
            MeshObjWriter::write(*_mesh, tetSurfaceMeshFile.c_str()); 

#ifndef USE_ADF
        _signedDistanceField.reset(
                DistanceFieldBuilder::BuildSignedClosestPointField(
                    tetSurfaceMeshFile.c_str(), 
                    _signedDistanceFieldResolution, 
                    _signedDistanceFieldFilePrefix.c_str()
                    )
                );
#else
        _signedDistanceField.reset(
                DistanceFieldBuilder::BuildAdaptiveDistanceField(
                    tetSurfaceMeshFile.c_str(),
                    _signedDistanceFieldFilePrefix.c_str(),
                    ADF_SUBDIVIDE_RADIUS,
                    ADF_MAX_OCTREE_LEVELS,
                    ADF_ERROR_TOLERANCE
                    )
                );
#endif
    }
    else // build from surface mesh
    {
        // first established paths for tet mesh, surface mesh, and correponding sdf
        const std::string meshFile = IO::AssembleFilePath(_workingDirectory, _objectPrefix) + ".obj"; 
        _signedDistanceFieldFilePrefix = meshFile + "." + std::to_string(_signedDistanceFieldResolution) + ".dist";

        if (!IO::ExistFile(meshFile))
            throw std::runtime_error("**ERROR** Surface mesh file not exist: " + meshFile); 
        
        _mesh = std::make_shared<TriangleMeshKDTree<REAL> >(); 
        //_mesh.reset(new TriangleMesh<REAL>());
        if (MeshObjReader::read(meshFile.c_str(), *_mesh, false, false, _meshScale)==SUCC_RETURN)
        {
            _mesh->generate_normals(); 
            _mesh->generate_mean_curvatures(); 
        }
        else
            throw std::runtime_error("**ERROR** Cannot read mesh from" + meshFile);

#ifndef USE_ADF
        _signedDistanceField.reset(
                DistanceFieldBuilder::BuildSignedClosestPointField(
                    meshFile.c_str(), 
                    _signedDistanceFieldResolution, 
                    _signedDistanceFieldFilePrefix.c_str()
                    )
                );
#else
        _signedDistanceField.reset(
                DistanceFieldBuilder::BuildAdaptiveDistanceField(
                    meshFile.c_str(),
                    _signedDistanceFieldFilePrefix.c_str(),
                    ADF_SUBDIVIDE_RADIUS,
                    ADF_MAX_OCTREE_LEVELS,
                    ADF_ERROR_TOLERANCE
                    )
                ); 
#endif
    }

    // compute mesh centroid in object space and cache it 
    _meshObjectCentroid = _mesh->ComputeCentroid(); 

    std::cout << "Read in TriangleMesh: \n"
              << " name:       " << GetMeshName() << "\n"
              << " #vertices:  " << _mesh->num_vertices() << "\n"
              << " #triangles: " << _mesh->num_triangles() << "\n";


    // build kd-tree for query
    std::dynamic_pointer_cast<TriangleMeshKDTree<REAL> >(_mesh)->BuildKDTree(); 
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
void FDTD_RigidObject::
ResetUnionBox()
{
    _bboxWorldUnion2Steps = _bboxWorld;
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
    position = _modelingTransformInverse*position.eval();
    return _signedDistanceField->distance(Vector3d(position[0],position[1],position[2])); 
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
    position = _modelingTransformInverse*position.eval();
    queriedNormal = _signedDistanceField->gradient(Conversions::ToVector3<double>(position));
    Eigen::Vector3d normal = Conversions::ToEigen<double>(queriedNormal); 
    normal = _modelingTransform.linear()*normal.eval(); 
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
    REAL bcValue = 0.0; 
    const SourceIterator sourceEnd = _vibrationalSources.end(); 
    for (SourceIterator it=_vibrationalSources.begin(); it!=sourceEnd; ++it) 
        bcValue += (*it)->Evaluate(boundaryPoint, boundaryNormal, time);

    return bcValue; 
}

//##############################################################################
// This function locates acceleration noise source and evaluate the analytical
// pressure for spheres.
//##############################################################################
REAL FDTD_RigidObject::
EvaluateAccelerationNoiseAnalytical(const Vector3d &listeningPoint, const REAL &time, const REAL &density, const REAL &soundSpeed, const REAL &sphereRadius)
{
    REAL bcValue = 0.0;
    const SourceIterator sourceEnd = _vibrationalSources.end(); 
    for (SourceIterator it=_vibrationalSources.begin(); it!=sourceEnd; ++it)
    {
        AccelerationNoiseVibrationalSource *anSource = dynamic_cast<AccelerationNoiseVibrationalSource*>((*it).get()); 
        if (anSource)
        {
            bcValue += anSource->EvaluatePressureAnalytical(listeningPoint, Vector3d(0,0,0), time, density, soundSpeed, sphereRadius); 
        }
        else 
        {
        }
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
// If reflection fails to push the point out of the obundary, then we fall back 
// to the nearest neighbour triangle search in order to establish a valid boundary 
// point and reflection. The checked condition is:
//  1. If image point is still inside the geometry. 
//
// If the original query point is already outside the boundary, this function
// will extend it in the normal direction and still perform a "reflection". 
// See below diagram.
//
//           /
// RP  OP  BP
//  o---o---o  inside boundary
//          |  
//          |
//##############################################################################
bool FDTD_RigidObject::
ReflectAgainstBoundary(const Vector3d &originalPoint, Vector3d &reflectedPoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &distanceTravelled)
{
    assert(_signedDistanceField!=nullptr); //&& (DistanceToMesh(originalPoint.x,originalPoint.y,originalPoint.z)<DISTANCE_TOLERANCE));

#if 0 // use sdf for normal query
    // find boundary point, normal at query point, and reflection point.
    NormalToMesh(originalPoint.x, originalPoint.y, originalPoint.z, erectedNormal);
    erectedNormal.normalize(); 
    distanceTravelled = DistanceToMesh(originalPoint.x, originalPoint.y, originalPoint.z);
    if (distanceTravelled > DISTANCE_TOLERANCE) // located outside the boundary already, push it further
    {
        boundaryPoint = originalPoint - erectedNormal * (distanceTravelled);
        reflectedPoint= originalPoint + erectedNormal * (distanceTravelled); 
    }
    else // inside the boundary, follow a similar procedure as the ghost cell method 
    {
        boundaryPoint = originalPoint - erectedNormal * (distanceTravelled);
        reflectedPoint= originalPoint - erectedNormal * (2.0*distanceTravelled); // want it to be outside the boundary
    }

#else // use kd-tree for normal query

    const Vector3d originalPointObject = WorldToObjectPoint(originalPoint); 
    int closestTriangleIndex; 
    Vector3d projectedPoint; // object space
    distanceTravelled = _mesh->ComputeClosestPointOnMesh(originalPointObject, boundaryPoint, closestTriangleIndex, projectedPoint); 
    boundaryPoint = ObjectToWorldPoint(boundaryPoint); 

    if (distanceTravelled < KD_NEAREST_TOLERANCE) // dont trust the result if lower than tolerance
    {
        // get closest triangle normal and push manually
        Vector3d t_normal = _mesh->triangle_normal(closestTriangleIndex); 
        t_normal.normalize(); 
        reflectedPoint = originalPointObject + t_normal * TRI_NORMAL_PUSH_DIST; 
        distanceTravelled = (TRI_NORMAL_PUSH_DIST/2.0); 

        // transform
        reflectedPoint = ObjectToWorldPoint(reflectedPoint); 
        erectedNormal = ObjectToWorldVector(t_normal); 
        erectedNormal.normalize(); 
        boundaryPoint = (originalPoint + reflectedPoint)/2.0;
    }
    else
    {
        const bool insideBoundary = DistanceToMesh(originalPoint.x, originalPoint.y, originalPoint.z) < 0 ? true : false;
        if (insideBoundary)
        {
            erectedNormal = boundaryPoint - originalPoint; // world space
            reflectedPoint = boundaryPoint + erectedNormal; 
        }
        else 
        {
            erectedNormal =-boundaryPoint + originalPoint; // world space
            reflectedPoint = originalPoint + erectedNormal; 
        }
        erectedNormal.normalize();  // keep the behavior same as the distance field based query
    }

#endif // if 0

#if DEBUG_WRITE_REFLECTION_ARROWS_INTERVAL > 0
#ifdef USE_OPENMP
#pragma omp critical
#endif
    {
        //const Vector3d &vertex = closestPoint; 
        //const Eigen::Vector3d vertexWorld_e = _modelingTransform * Eigen::Vector3d(vertex.x, vertex.y, vertex.z); 
        //const Vector3d vertexWorld(vertexWorld_e[0], vertexWorld_e[1], vertexWorld_e[2]); 
        // write these special points for debugging purpose 
        _debugArrowStart.push_back(originalPoint); 
        _debugArrowNormal.push_back(reflectedPoint - originalPoint); 
    }
#endif

    const REAL newDistance = DistanceToMesh(reflectedPoint);
    const bool reflectSuccess = (newDistance >= DISTANCE_TOLERANCE); 
    return reflectSuccess;
}

//##############################################################################
// This function is now deprecated. call ReflectAgainstBoundary directly.
//##############################################################################
bool FDTD_RigidObject::
FindImageFreshCell(const Vector3d &currentPoint, Vector3d &imagePoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &distanceTravelled)
{
    return ReflectAgainstBoundary(currentPoint, imagePoint, boundaryPoint, erectedNormal, distanceTravelled); 
}

//##############################################################################
//##############################################################################
Vector3d FDTD_RigidObject::
MeshCentroid()
{
    return ObjectToWorldPoint(_meshObjectCentroid); 
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

//##############################################################################
//##############################################################################
void FDTD_RigidObject::
WriteDebugArrow(const std::string &file)
{
    std::ofstream of(file.c_str()); 
    for (size_t idx=0; idx<_debugArrowStart.size(); ++idx)
    {
        of << _debugArrowStart.at(idx).x << " " << _debugArrowStart.at(idx).y << " " << _debugArrowStart.at(idx).z << " " 
           << _debugArrowNormal.at(idx).x << " " << _debugArrowNormal.at(idx).y << " " << _debugArrowNormal.at(idx).z << std::endl;
    }
    of.close();
    _debugArrowStart.clear(); 
    _debugArrowNormal.clear(); 
}

//##############################################################################
//##############################################################################
void FDTD_RigidObject::
ClearDebugArrow()
{
    _debugArrowStart.clear(); 
    _debugArrowNormal.clear(); 
}
