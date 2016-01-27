#ifndef TETRAHEDRALMESH_H 
#define TETRAHEDRALMESH_H 

#include <Eigen/Dense> 
#include <iostream>
#include "KDTree.hpp" 

class VecDistSqr
{
public:
    template<typename MatrixType>
    double operator() (const MatrixType &v1, const MatrixType &v2) const
    {
        return (v1-v2).squaredNorm();
    }
};

class VecAcc
{
public:
    typedef double result_type;

    template<typename MatrixType>
    result_type operator() (const MatrixType &v1, const size_t N) const
    {
        return v1(N);
    }
};

/* 
 * stroe some extra data on the tetrahedron
 */
struct TetrahedronExtraData
{
    int openfoamCellIndex; // memorize the openfoam cell index this tetrahedron is associated to
}; 


/* 
 * tetrahedral primitive
 */
class Tetrahedron
{
    private: 

        Eigen::Vector3d _p0; 
        Eigen::Vector3d _p1; 
        Eigen::Vector3d _p2; 
        Eigen::Vector3d _p3; 

        TetrahedronExtraData    _extraData;

    public: 

        inline Eigen::Vector3d Centroid() const { return (_p0+_p1+_p2+_p3)/4.0; }  
        inline void SetExtraData(const TetrahedronExtraData &extraData) { _extraData = extraData; }
        inline TetrahedronExtraData& GetExtraData() { return _extraData; }

        bool Inside(const Eigen::Vector3d &pos) const; 
        bool Inside(const Eigen::Vector3d &pos, Eigen::Vector4d &weights) const; 

        Tetrahedron(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, const Eigen::Vector3d &p3) : 
            _p0(p0), 
            _p1(p1), 
            _p2(p2), 
            _p3(p3) { 
            } 

        friend std::ostream& operator << (std::ostream &os, const Tetrahedron &t) 
        {
            os << "Tet: \n" 
               << t._p0.transpose() << std::endl
               << t._p1.transpose() << std::endl
               << t._p2.transpose() << std::endl
               << t._p3.transpose();

            return os; 
        }
};

/* 
 * tetrahedral mesh 
 */
class TetrahedronMesh
{
    public: 
        typedef KDTree<3, Eigen::Vector3d, VecDistSqr, VecAcc, Eigen::aligned_allocator<Eigen::Vector3d>> KdTree; 
        typedef std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d>> CentroidVector; 

    private: 
        KdTree                      _tetrahedronTree; 
        CentroidVector              _tetrahedronsCentroids; 
        std::vector<Tetrahedron>    _tetrahedrons; 

    public: 

        inline int N_tetrahedrons() { return _tetrahedrons.size(); }
        inline Tetrahedron& Get(const int &ind) { return _tetrahedrons[ind]; }
        int FindNearestTetrahedron(const Eigen::Vector3d &position); 
        int FindEnclosedTetrahedron(const Eigen::Vector3d &position, const int &N_neighbors=100); 
        void AddTetrahedron(const Tetrahedron &tet);
        void BuildTree(); 

        TetrahedronMesh(){ };

}; 



#endif 
