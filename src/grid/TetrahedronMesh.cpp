#include "TetrahedronMesh.h" 


bool Tetrahedron::Inside(const Eigen::Vector3d &p) const 
{
    Eigen::Vector4d weights; 
    const bool result = Inside(p,weights); 
    return result;
}

bool Tetrahedron::Inside(const Eigen::Vector3d &p, Eigen::Vector4d &weights) const
{

    Eigen::Matrix3d A; 
    A.col(0) = _p0 - _p3; 
    A.col(1) = _p1 - _p3; 
    A.col(2) = _p2 - _p3; 

    Eigen::Vector3d b = p - _p3; 

    const double detA = A.determinant(); 
    double sumLambdas = 0; 

    for (int ii=0; ii<3; ii++) 
    { 
        Eigen::Matrix3d D = A; 
        D.col(ii) = b; 
        const double detD = D.determinant(); 
        weights(ii) = detD / detA; 
        if (weights(ii)<0.0 || weights(ii)>1.0) // not within tet
            return false; 
        sumLambdas += weights(ii); 
    } 

    weights(3) = 1.0 - sumLambdas;

    if (weights(3)<0.0|| weights(3)>1.0)
        return false; 

    return true;
}

void TetrahedronMesh::AddTetrahedron(const Tetrahedron &tet) 
{
    _tetrahedrons.push_back(tet); 
    _tetrahedronsCentroids.push_back(tet.Centroid());
}

void TetrahedronMesh::BuildTree()
{
    std::cout << "building kd-tree for centroid nearest neighbour search ..." << std::flush; 
    if (_tetrahedrons.size()>0) 
        _tetrahedronTree.reinitialize( _tetrahedronsCentroids.data(), _tetrahedronsCentroids.size() );
    std::cout << " OK" << std::endl;
}

int TetrahedronMesh::FindNearestTetrahedron(const Eigen::Vector3d &position)
{
    return _tetrahedronTree.find_nearest(position);
}

int TetrahedronMesh::FindEnclosedTetrahedron(const Eigen::Vector3d &position, const int &N_neighbors) 
{
    std::vector<int> indicies;
    _tetrahedronTree.find_nearest(position, N_neighbors, indicies); 
    //_tetrahedronTree.find_in_ball(position, 0.01, indicies); 

    int enclosedTetrahedronIndex = -1; 
    for (size_t ii=0; ii<indicies.size(); ii++) 
    {
        Eigen::Vector4d weight;
        if (_tetrahedrons[indicies[ii]].Inside(position, weight))
            enclosedTetrahedronIndex = indicies[ii]; 
    }

    return enclosedTetrahedronIndex; 
}
