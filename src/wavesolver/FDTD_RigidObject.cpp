#include <wavesolver/FDTD_RigidObject.h> 
#include <io/TglMeshReader.hpp>
#include <utils/SimpleTimer.h>

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

    return _signedDistanceField->distance(position); 
}

//##############################################################################
//##############################################################################
void FDTD_RigidObject::
NormalToMesh(const double &x, const double &y, const double &z, Vector3d &queriedNormal)
{
    if (!_signedDistanceField)
        throw std::runtime_error("**ERROR** distance field not built.");

    if (!_bboxWorld.Inside(x,y,z,1.1))
        return std::numeric_limits<REAL>::lowest();

    Eigen::Vector3d position(x,y,z); 
    position = _modelingTransformInverse*position;

    queriedNormal = _signedDistanceField->gradient(position);
}

//##############################################################################
//##############################################################################
void FDTD_RigidObject::
PrintBoundingBox()
{
    std::cout << "minBound = " << _bboxWorld.minBound << std::endl; 
    std::cout << "maxBound = " << _bboxWorld.maxBound << std::endl; 
}
