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
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrixXd; 

    protected:
        std::shared_ptr<VlKDForest>         _nnForest; 
        std::vector<Vector3<T>>             _triangleCentroids; 

    public: 
        ~TriangleMeshKDTree()
        {
            //vl_kdforest_delete(_nnForest.get());
        }

        virtual T FindKNearestTriangles(const int &k, const Vector3<T> &point, std::vector<int> &triangleIndices) const; 
        virtual T FindNearestTriangle(const Vector3<T> &point, int &triangleIndex) const; 
        virtual int FindTrianglesWithinBall(const Vector3<T> &point, const T &radius, std::set<int> &triangleIndices) const; 
        inline const Vector3<T> &TriangleCentroid(const int &t_idx) const {return _triangleCentroids.at(t_idx);}
        void BuildKDTree();

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
    //_triangleCentroids.resize(N_triangles, 3); 
    _triangleCentroids.resize(N_triangles); 
    for (int t_idx=0; t_idx<N_triangles; ++t_idx)
    {
        const Tuple3ui &indices  = this->m_triangles.at(t_idx); 
        const Point3<T> &vertex0 = this->m_vertices.at(indices.x); 
        const Point3<T> &vertex1 = this->m_vertices.at(indices.y); 
        const Point3<T> &vertex2 = this->m_vertices.at(indices.z); 
        const Point3<T> centroid = (vertex0 + vertex1 + vertex2)/3.0;
        _triangleCentroids.at(t_idx) = centroid;
    }

    // build forest
    _nnForest.reset(vl_kdforest_new(VL_TYPE_DOUBLE, 3, 1, VlDistanceL2)); 
    vl_kdforest_build(_nnForest.get(), N_triangles, &(_triangleCentroids[0])); 
}

//##############################################################################
// @param k Number of nearest neighbours
// @param point Sample point position
// @param triangleIndices Output triangles indices
// @return smallest distance
//##############################################################################
template <typename T> 
T TriangleMeshKDTree<T>::
FindKNearestTriangles(const int &k, const Vector3<T> &point, std::vector<int> &triangleIndices) const
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
// This function is a special case of finding k nearest neighbours. However for
// performace I copied some of the codes.
// @param k Number of nearest neighbours
// @param point Sample point position
// @param triangleIndices Output triangles indices
// @return smallest distance
//##############################################################################
template <typename T> 
T TriangleMeshKDTree<T>::
FindNearestTriangle(const Vector3<T> &point, int &triangleIndex) const
{
#if defined(USE_BOOST) && defined(DEBUG_KDTREE)
    boost::timer::auto_cpu_timer timer("Boost timer: KDTree Query takes %w seconds\n" );
#endif
    VlKDForestNeighbor neighbours; 
    vl_kdforest_query(_nnForest.get(), &neighbours, 1, &point); 
    triangleIndex = neighbours.index; 
    return neighbours.distance; 
}

//##############################################################################
// Function FindTrianglesWithinBallNaive
//   @param point Sample point position
//   @param radius ball radius 
//   @param triangleIndices output field that stores all triangles within the ball
//   @return number of triangles
//##############################################################################
template <typename T> 
int TriangleMeshKDTree<T>::
FindTrianglesWithinBall(const Vector3<T> &point, const T &radius, std::set<int> &triangleIndices) const
{
    triangleIndices.clear(); 
    if (radius<=0)
        return 0; 
#if 0 // naive
    const int N = this->num_triangles(); 
    for (int ii=0; ii<N; ++ii)
        if ((this->TriangleCentroid(ii)-point).length() <= radius)
            triangleIndices.insert(ii); 
#else
    const int increment = 200; 
    int knn = increment; 
    std::vector<int> indices; 
    typedef std::vector<int>::iterator It; 
    while (true) 
    {
        indices.clear(); 
        FindKNearestTriangles(knn, point, indices);
        const int farest = std::min<int>(indices.size(), this->num_triangles()); 
        if (indices.size()>=this->num_triangles() || 
            (this->TriangleCentroid(indices.at(farest-1)) - point).length() > radius)
        {
            It it = indices.begin();
            for (; it!=indices.end(); ++it)
                if ((this->TriangleCentroid(*it)-point).length() > radius)
                    break; 
            if (indices.size()>=this->num_triangles())
                indices.erase(indices.begin()+this->num_triangles(), indices.end()); 
            else
                indices.erase(it, indices.end()); 
            triangleIndices = std::set<int>(indices.begin(), indices.end()); 
            break; 
        }
        knn += increment;
    }
#endif
    return triangleIndices.size();
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
    const T boundingRadius = this->boundingSphereRadius(Point3<T>(0, 0, 0)); 
    const int N_triangles = this->triangles().size(); 
    srand(time(NULL)); 
    std::vector<int> triangleIndices;

    int misClassifyIndices = 0; 
    T misClassifyError = 0.0;
    for (int s_idx=0; s_idx<N_samples; ++s_idx)
    {
        const T x = static_cast<T>(rand())/static_cast<T>(RAND_MAX) * boundingRadius;
        const T y = static_cast<T>(rand())/static_cast<T>(RAND_MAX) * boundingRadius;
        const T z = static_cast<T>(rand())/static_cast<T>(RAND_MAX) * boundingRadius;
        const Vector3<T> sample(x, y, z); 

        // do a brute force search to find reference;
        T minDistance = std::numeric_limits<T>::max(); 
        int index = -1; 
        for (int t_idx=0; t_idx<N_triangles; ++t_idx)
        {
            const T distance = (TriangleCentroid(t_idx) - sample).lengthSqr(); 
            if (distance < minDistance) 
            {
                minDistance = distance; 
                index = t_idx;
            }
        }

        // knn query
        const T closestDistance = FindKNearestTriangles(k, sample, triangleIndices);

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
