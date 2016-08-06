#ifndef TRIANGLE_MESH_KD_TREE_H
#define TRIANGLE_MESH_KD_TREE_H 

#include <linearalgebra/Vector3.hpp>
#include <geometry/TriangleMesh.hpp> 
#include <macros.h>
// VLFeat is C library, so need extern
extern "C" {
#include <generic.h>
#include <kdtree.h>
}

#ifdef USE_BOOST
#include <boost/timer/timer.hpp>
#endif

// uncomment if want the timer for each query
//#define DEBUG_KDTREE

//##############################################################################
// This class will put the triangle mesh class into a kd-tree and support
// subsequent query about nearest neighbours. 
//
// Currently, the kd-tree building and querying is done via external library
// VLFeat (http://www.vlfeat.org/index.html). Documentation for the kd-tree
// operations can be found here: http://www.vlfeat.org/api/kdtree.html 
//##############################################################################
template <typename T> 
class TriangleMeshKDTree : public TriangleMesh<T>
{
    public: 
        typedef Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrixXd; 

    protected:
        std::shared_ptr<VlKDForest>         _nnForest; 
        RowMajorMatrixXd                    _triangleCentroids; 

    public: 
        ~TriangleMeshKDTree()
        {
            //vl_kdforest_delete(_nnForest.get());
        }

        void BuildKDTree();
        REAL FindKNearestTriangles(const int &k, const Vector3d &point, std::vector<int> &triangleIndices); 

        ///// debugging methods /////
        void TestKDTree(const int &k); 
};

//##############################################################################
//##############################################################################
template <typename T> 
void TriangleMeshKDTree<T>::
BuildKDTree()
{
#ifdef USE_BOOST
    boost::timer::auto_cpu_timer timer("Boost timer: Building KDTree for mesh takes %w seconds\n" );
#endif
    if (_nnForest)
        return; 

    // first compute centroids for each triangles and stored them
    const int N_triangles = this->m_triangles.size(); 
    _triangleCentroids.resize(N_triangles, 3); 
    for (int t_idx=0; t_idx<N_triangles; ++t_idx)
    {
        const Tuple3ui &indices  = this->m_triangles.at(t_idx); 
        const Point3<T> &vertex0 = this->m_vertices.at(indices.x); 
        const Point3<T> &vertex1 = this->m_vertices.at(indices.y); 
        const Point3<T> &vertex2 = this->m_vertices.at(indices.z); 
        const Point3<T> centroid = (vertex0 + vertex1 + vertex2)/3.0;
        _triangleCentroids(t_idx, 0) = centroid.x;
        _triangleCentroids(t_idx, 1) = centroid.y;
        _triangleCentroids(t_idx, 2) = centroid.z;
    }

    // build forest
    _nnForest.reset(vl_kdforest_new(VL_TYPE_DOUBLE, 3, 1, VlDistanceL2)); 
    vl_kdforest_build(_nnForest.get(), N_triangles, _triangleCentroids.data()); 
}

//##############################################################################
// @param k Number of nearest neighbours
// @param point Sample point position
// @param triangleIndices Output triangles indices
// @return smallest distance
//##############################################################################
template <typename T> 
REAL TriangleMeshKDTree<T>::
FindKNearestTriangles(const int &k, const Vector3d &point, std::vector<int> &triangleIndices)
{
#if defined(USE_BOOST) && defined(DEBUG_KDTREE)
    boost::timer::auto_cpu_timer timer("Boost timer: KDTree Query takes %w seconds\n" );
#endif
    VlKDForestNeighbor neighbours[k]; 
    vl_kdforest_query(_nnForest.get(), neighbours, k, &point); 
    triangleIndices.resize(k); 
    for (int nn_idx=0; nn_idx<k; ++nn_idx)
        triangleIndices.at(nn_idx) = neighbours[nn_idx].index; 
    return neighbours[0].distance; 
}

//##############################################################################
//##############################################################################
template <typename T> 
void TriangleMeshKDTree<T>::
TestKDTree(const int &k)
{
#ifdef USE_BOOST
    boost::timer::auto_cpu_timer timer("Boost timer: Testing KDTree takes %w seconds\n" );
#endif
    const int N_samples = 200; // test on 200 queries
    const REAL boundingRadius = this->boundingSphereRadius(Point3<T>(0, 0, 0)); 
    const int N_triangles = this->triangles().size(); 
    srand(time(NULL)); 
    std::vector<int> triangleIndices;

    int misClassifyIndices = 0; 
    REAL misClassifyError = 0.0;
    for (int s_idx=0; s_idx<N_samples; ++s_idx)
    {
        const REAL x = static_cast<REAL>(rand())/static_cast<REAL>(RAND_MAX) * boundingRadius;
        const REAL y = static_cast<REAL>(rand())/static_cast<REAL>(RAND_MAX) * boundingRadius;
        const REAL z = static_cast<REAL>(rand())/static_cast<REAL>(RAND_MAX) * boundingRadius;
        const Vector3d sample(x, y, z); 
        RowMajorMatrixXd sample_e(1, 3); 
        sample_e << x, y, z; 

        // do a brute force search to find reference;
        REAL minDistance = std::numeric_limits<REAL>::max(); 
        int index = -1; 
        for (int t_idx=0; t_idx<N_triangles; ++t_idx)
        {
            const REAL distance = (_triangleCentroids.row(t_idx) - sample_e).squaredNorm(); 
            if (distance < minDistance) 
            {
                minDistance = distance; 
                index = t_idx;
            }
        }

        // knn query
        const REAL closestDistance = FindKNearestTriangles(k, sample, triangleIndices);

        // compute and store testing information
        if (index != triangleIndices[0] || !EQUAL_FLOATS(minDistance, closestDistance))
        {
            misClassifyIndices += 1;
            misClassifyError = max(misClassifyError, abs(minDistance-closestDistance)); 
        }
    }
    std::cout << "====================================================\n";
    std::cout << "KD Tree test results: \n"; 
    std::cout << " Queried number of nearest neighbours (exact) = " << k << "\n";
    std::cout << " Mis-classified triangle indices = " << misClassifyIndices << "\n";
    std::cout << " Mis-classified distance error = "   << misClassifyError << "\n";
    if (misClassifyIndices == 0)
        std::cout << "**TEST PASSED\n";
    else
        std::cout << "**TEST FAILED\n";
    std::cout << "====================================================";
    std::cout << std::endl;
}

#endif
