#include <wavesolver/MAC_Grid.h>
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
Initialize(const bool &buildFromTetMesh, bool needsCurvature)
{
    assert(_parsed); 

    _disableEvals = false;
    if (_signedDistanceFieldResolution > 0)
    {
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
                _mesh = std::make_shared<TriangleMeshGraph<REAL> >(); 
                //_mesh.reset(new TriangleMesh<REAL>()); 
                tetMesh->extract_surface(_mesh.get()); 
                _mesh->generate_normals(); 
                _mesh->update_vertex_areas(); 
                _tetMeshIndexToSurfaceMesh = std::make_shared<TetMeshIndexToSurfaceMesh>(); 
                _tetMeshIndexToSurfaceMesh->ReadFromGeoFile(geoFile); 
                if (_tetMeshIndexToSurfaceMesh->N_surfaceVertices() != _mesh->num_vertices())
                    throw std::runtime_error("**ERROR** geo file has different number of surface vertices than the surface mesh from tet mesh");
                else 
                    std::cout << " Surface mesh and Tet-Surface mapping built.\n";
                _volume = tetMesh->total_volume(); // cache volume
                tetMesh->inertia_tensor(_volumeInertiaTensor, _volumeCenter); 
                _hasVolume = true;
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
            
            _mesh = std::make_shared<TriangleMeshGraph<REAL> >(); 
            //_mesh.reset(new TriangleMesh<REAL>());
            if (MeshObjReader::read(meshFile.c_str(), *_mesh, false, false, _meshScale)==SUCC_RETURN)
            {
                _mesh->generate_normals(); 
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
    }
    else  // sdf resolution <= 0
    {
        // first established paths for tet mesh, surface mesh, and correponding sdf
        const std::string meshFile = IO::AssembleFilePath(_workingDirectory, _objectPrefix) + ".obj"; 
        std::cout << "FIXME debug " << meshFile << std::endl;
        if (!IO::ExistFile(meshFile))
            throw std::runtime_error("**ERROR** Surface mesh file not exist: " + meshFile); 
        _mesh = std::make_shared<TriangleMeshGraph<REAL> >(); 
        if (MeshObjReader::read(meshFile.c_str(), *_mesh, false, false, _meshScale)==SUCC_RETURN)
        {
            _mesh->generate_normals(); 
        }
        else
            throw std::runtime_error("**ERROR** Cannot read mesh from" + meshFile);
    
    }

    // compute mesh centroid in object space and cache it 
    _meshObjectCentroid = _mesh->ComputeCentroid(); 

    // get curvatures info
    const int N_vertices = _mesh->vertices().size(); 
    REAL maxCurvature = std::numeric_limits<REAL>::min(); 
    REAL minCurvature = std::numeric_limits<REAL>::max(); 

    if (needsCurvature)
    {
        _mesh->generate_mean_curvatures(); 
        const std::vector<REAL> *meanCurvatures = _mesh->mean_curvatures();
        for (int v_idx=0; v_idx<N_vertices; ++v_idx)
        {
            maxCurvature = std::max<REAL>(maxCurvature, meanCurvatures->at(v_idx)); 
            minCurvature = std::min<REAL>(minCurvature, meanCurvatures->at(v_idx)); 
        }
    }

    std::cout << "Read in TriangleMesh: \n"
              << " Name           : " << GetMeshName() << "\n"
              << " Num vertices   : " << _mesh->num_vertices() << "\n"
              << " Num triangles  : " << _mesh->num_triangles() << "\n"
              << " Curvature range: [" << minCurvature << ", " << maxCurvature  << "] " << std::endl;

    // build kd-tree and graph for query
    std::dynamic_pointer_cast<TriangleMeshKDTree<REAL> >(_mesh)->BuildKDTree(); 
    _meshGraph = std::dynamic_pointer_cast<TriangleMeshGraph<REAL> >(_mesh); 
    assert(_solverSettings);
    int largestTriangle; 
    const REAL d_max = _solverSettings->cellSize; 
    const REAL t_max = sqrt(_mesh->largest_triangle_area(largestTriangle)); 
    //const REAL d_max=0.0;
    //const REAL t_max=0.0;
    const double knnradius = d_max + t_max; 
    const std::string meshGraphFile = IO::AssembleFilePath(_workingDirectory, _objectPrefix) + 
        "." + std::to_string(knnradius) + ".meshgraph"; 
    _meshGraph->BuildGraph(meshGraphFile, knnradius); 
    UpdateBoundingBox();
}

//##############################################################################
//##############################################################################
void FDTD_RigidObject::
Reinitialize(const std::string& objPrefix, const bool buildFromTetMesh, bool needsCurvature)
{
    _objectPrefix = objPrefix;
    Initialize(buildFromTetMesh, needsCurvature);
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
        Vector3<REAL> minBound( D_INF,  D_INF,  D_INF); 
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
        if( minBound.x == D_INF )
            throw std::runtime_error("**ERROR** bounding box not set in UpdateBoundingBox");
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
ApplyScale(const REAL scale)
{
    FDTD_MovableObject::ApplyScale(scale); // scale the transformation
    _meshScale *= scale; 
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
    if (!_bboxWorld.Inside(x, y, z, AABB_CHECK_TOLERANCE_SCALE))
        return std::numeric_limits<REAL>::max();
    Eigen::Vector3d position(x,y,z); 
    position = _modelingTransformInverse*position.eval();
    REAL d; 
    if (!_signedDistanceField)
    {
        return std::numeric_limits<REAL>::max(); 
    }
    else
    {
        d = _signedDistanceField->distance(
              Vector3d(position[0],
                       position[1],
                       position[2])); 
    }

    return d * _meshScale; 
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
    if (!_bboxWorld.Inside(x,y,z,1.1))
    {
        const REAL limit = std::numeric_limits<REAL>::max(); 
        queriedNormal = Vector3d(limit, limit, limit); 
        return false;
    }

    Eigen::Vector3d position(x,y,z); 
    position = _modelingTransformInverse*position.eval();
    if (!_signedDistanceField)
    {
        return std::numeric_limits<REAL>::max(); 
    }
    else
    {
        queriedNormal = _signedDistanceField->gradient(Conversions::ToVector3<double>(position));
        Eigen::Vector3d normal = Conversions::ToEigen<double>(queriedNormal); 
        normal = _modelingTransform.linear()*normal.eval(); 
        queriedNormal = Conversions::ToVector3(normal); 
    }
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
EvaluateBoundaryAcceleration(const Vector3d &boundaryPoint, const Vector3d &boundaryNormal, const REAL &time, const int &hintTriangle)
{
    REAL bcValue = 0.0; 
    if (_disableEvals) return 0.0;
    const SourceIterator sourceEnd = _vibrationalSources.end(); 
    for (SourceIterator it=_vibrationalSources.begin(); it!=sourceEnd; ++it) 
        bcValue += (*it)->Evaluate(boundaryPoint, boundaryNormal, time, hintTriangle);

    return bcValue; 
}

//##############################################################################
//##############################################################################
REAL FDTD_RigidObject::
EvaluateBoundaryAcceleration(const int &vertexID, const Vector3d &vertexNormal, const REAL &time)
{
    REAL bcValue = 0.0; 
    if (_disableEvals) return 0.0;
    const SourceIterator sourceEnd = _vibrationalSources.end(); 
    for (SourceIterator it=_vibrationalSources.begin(); it!=sourceEnd; ++it) 
        bcValue += (*it)->Evaluate(vertexID, vertexNormal, time);

    return bcValue; 
}

//##############################################################################
//##############################################################################
Vector3d FDTD_RigidObject::
EvaluateBoundaryAcceleration(const int &vertexID, const REAL &time)
{
    Vector3d bcValue; 
    if (_disableEvals) return bcValue;
    const SourceIterator sourceEnd = _vibrationalSources.end(); 
    for (SourceIterator it=_vibrationalSources.begin(); it!=sourceEnd; ++it) 
        bcValue += (*it)->Evaluate(vertexID, time);

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
    if (_disableEvals) return 0.0;
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
    if (_disableEvals) return 0.0;
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
int FDTD_RigidObject::
ReflectAgainstBoundary(const Vector3d &originalPoint, Vector3d &reflectedPoint, Vector3d &boundaryPoint, Vector3d &erectedNormal, REAL &distanceTravelled, const int &startFromTriangle)
{
#if 0 // use sdf for normal query
    assert(_signedDistanceField); //&& (DistanceToMesh(originalPoint.x,originalPoint.y,originalPoint.z)<DISTANCE_TOLERANCE));
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
        distanceTravelled = -distanceTravelled; // make it positive
    }
    return -1; 

#else // use kd-tree for normal query

    const Vector3d originalPointObject = WorldToObjectPoint(originalPoint); 
    int closestTriangleIndex; 
    Vector3d projectedPoint; // object space
    if (startFromTriangle<0) 
    {
        distanceTravelled = _mesh->ComputeClosestPointOnMesh(originalPointObject, boundaryPoint, closestTriangleIndex, projectedPoint, 100); 
    }
    else
    {
        distanceTravelled = _meshGraph->ComputeClosestPointOnMesh(startFromTriangle, originalPointObject, boundaryPoint, closestTriangleIndex, projectedPoint, 0.99, 100); 
        // uncomment if want to test between graph-search and kdtree-search results
        //int graph = closestTriangleIndex; 
        //const Vector3d bpGraph = boundaryPoint; 
        //distanceTravelled = _mesh->ComputeClosestPointOnMesh(originalPointObject, boundaryPoint, closestTriangleIndex, projectedPoint); 
        //if (graph != closestTriangleIndex) 
        //{
        //    std::cerr << "DIFFERENCE IN SEARCH: " << graph << " <-> " << closestTriangleIndex << std::endl; 
        //    std::cerr << " point       = " << originalPoint << std::endl; 
        //    std::cerr << " bp_graph    = " << ObjectToWorldPoint(bpGraph)       << std::endl; 
        //    std::cerr << " bp_kdtre    = " << ObjectToWorldPoint(boundaryPoint) << std::endl; 
        //}
    }
    boundaryPoint = ObjectToWorldPoint(boundaryPoint); 
    if (fabs(distanceTravelled) < KD_NEAREST_TOLERANCE) // dont trust the result if lower than tolerance
    {
        // get closest triangle normal and push manually
        Vector3d t_normal = _mesh->triangle_normal(closestTriangleIndex); 
        t_normal.normalize(); 
        reflectedPoint = originalPointObject + t_normal * 2.0*KD_NEAREST_TOLERANCE; 
        distanceTravelled = -KD_NEAREST_TOLERANCE; 

        // transform
        reflectedPoint = ObjectToWorldPoint(reflectedPoint); 
        erectedNormal = ObjectToWorldVector(t_normal); 
        erectedNormal.normalize(); 
        boundaryPoint = (originalPoint + reflectedPoint)/2.0;
    }
    else
    {
        const bool insideBoundary = (distanceTravelled < 0 ? true : false);
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
        erectedNormal.normalize();
        // check the normal compared to the triangle normal, they should always be pointing
        // in the same direction
        const Vector3d t_normal = ObjectToWorldVector(_mesh->triangle_normal(closestTriangleIndex)); 
        const REAL sgn = erectedNormal.dotProduct(t_normal)/t_normal.norm(); 
        if (sgn <= 0.0) 
        {
            throw std::runtime_error("**ERROR** triangle normal and erected normal are in the opposite direction."); 
        }
    }

#endif // if 0

#if DEBUG_WRITE_REFLECTION_ARROWS_INTERVAL > 0
#ifdef USE_OPENMP
#pragma omp critical
#endif
    {
        // write these special points for debugging purpose 
        _debugArrowStart.push_back(originalPoint); 
        _debugArrowNormal.push_back(reflectedPoint - originalPoint); 
    }
#endif
    return closestTriangleIndex;
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
void FDTD_RigidObject::
SetRigidBodyTransform(const Point3d &newCOM, const Quaternion<REAL> &quaternion)
{
    const Point3d &restCOM = _volumeCenter; 
    Vector3d rotationAxis; 
    const REAL rotationAngle = quaternion.toAxisRotR(rotationAxis);
    const Eigen::AngleAxisd rotation(rotationAngle, Eigen::Vector3d(rotationAxis.x, rotationAxis.y, rotationAxis.z));
    //const Eigen::Quaterniond rotation(quaternion.w, quaternion.v.x, quaternion.v.y, quaternion.v.x); 
    const Eigen::Vector3d restCOM_e = Eigen::Vector3d(restCOM.x, restCOM.y, restCOM.z); 
    const Eigen::Vector3d newCOM_e = Eigen::Vector3d(newCOM.x, newCOM.y, newCOM.z); 

    // reset transformation
    _modelingTransform.setIdentity();
    _modelingTransform.pretranslate(-restCOM_e);
    _modelingTransform.prerotate(rotation);
    _modelingTransform.pretranslate(newCOM_e); 
    _modelingTransformInverse = _modelingTransform.inverse(); 
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
REAL FDTD_RigidObject::
GetEarliestEventTime(const REAL &startTime) const
{
    REAL t = std::numeric_limits<REAL>::max(); 
    for (const auto &src : _vibrationalSources)
    {
        t = std::min(t, src->EarliestEventTime(startTime));
    }
    return t; 
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

//##############################################################################
//##############################################################################
int FDTD_RigidObject::
FindLowestVertex(const int &dimension, Vector3d &position)
{
    const std::vector<Point3<REAL> > &vertices = _mesh->vertices(); 
    const int N_vertices = vertices.size(); 
    int minIndex = -1; 
    REAL minValue = std::numeric_limits<REAL>::max(); 
    for (int v_idx=0; v_idx<N_vertices; ++v_idx)
    {
        const REAL value = ObjectToWorldPoint(vertices.at(v_idx))[dimension]; 
        if (value < minValue)
        {
            minIndex = v_idx; 
            minValue = value; 
        }
    }
    assert(minIndex != -1); 
    position = ObjectToWorldPoint(vertices.at(minIndex)); 
    std::cout << "vertex " << minIndex << " has the minimum position at dimension " << dimension << ": " << minValue << std::endl;
    return minIndex;
}
