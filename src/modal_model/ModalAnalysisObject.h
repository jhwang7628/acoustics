#ifndef MODAL_ANALYSIS_OBJECT_H 
#define MODAL_ANALYSIS_OBJECT_H 
#include <iostream>
#include <memory>
#include <deformable/ModeData.h> 
#include <modal_model/ModalODESolver.h> 
#include <modal_model/ModalMaterial.h>
#include <modal_model/SparseModalEncoder.h>
#include <geometry/TetMeshIndexToSurfaceMesh.h>
#include <geometry/TriangleMesh.hpp>
#include <Eigen/Dense>

//##############################################################################
// Class that handles the modal object
//
// Notations in the comments: 
//  N_V: number of volume tet mesh vertices
//  N_S: number of surface triangle mesh vertices
//  M: number of modes
//##############################################################################
class ModalAnalysisObject
{
    public: 
        typedef std::shared_ptr<ModalODESolver> ModalODESolverPtr; 
        typedef std::shared_ptr<TriangleMesh<REAL> > TriangleMeshPtr;
        enum IDType{VOL_UNCULL=0, SURF_CULLED=1}; // if read directly from *.modes its 0, after calling CullNonSurfaceModeShapes once should be set to 1

    private: 
        bool        _modeFileSet; 
        std::string _modeFile; 

        // volume <-> surface conversion
        std::shared_ptr<TetMeshIndexToSurfaceMesh> _tetMeshIndexToSurfaceMesh; 
        IDType                                     _idType; 

    protected: 
        Eigen::VectorXd _eigenValues; // eigenvalues produced by modal analysis (omega squared * density)
        Eigen::MatrixXd _eigenVectors; // dim: 3N_V * M before culling, 3N_S * M after culling
        Eigen::MatrixXd _eigenVectorsNormal; // eigenvectors projected to vertex normal direction. dim: N_S * M

        std::shared_ptr<ModalMaterial> _material;
        std::vector<ModalODESolverPtr> _modalODESolvers; 
        REAL _ODEStepSize; 
        REAL _time;  // corresponds to qNew

        SparseModalEncoderPtr _modalAccEncoder; 

    public: 
        ModalAnalysisObject()
            : _modeFileSet(false), _idType(VOL_UNCULL), _time(0.0)
        {} 
        ModalAnalysisObject(const std::string &modeFile) 
            : _modeFileSet(true), _modeFile(modeFile), _idType(VOL_UNCULL), _time(0.0)
        {}

        inline int N_Modes() const {return _eigenVectors.cols();}
        inline int N_DOF() const {return _eigenVectors.rows();}
        inline int N_vertices() const {return N_DOF()/3;}
        inline REAL GetODESolverTime(){return _time;}
        inline void SetModeFile(const std::string &modeFile){_modeFile = modeFile; _modeFileSet = true;}
        inline REAL GetModeFrequency(const int modeIndex) const 
        {
            return (modeIndex < _eigenValues.size() ? 
                    sqrt(_eigenValues(modeIndex) * _material->inverseDensity) / (2.0*M_PI) : 
                    0.0);
        }
        inline bool IDTypeIsSurf(){return _idType==SURF_CULLED;}
        inline std::shared_ptr<ModalMaterial> GetMaterial(){return _material;}
        void ReadModeFromFile(); 
        // Get U^T f, where U is the eigenvector (modal) matrix. Input vertexID should be in zero-based, tet mesh id. 
        //void GetForceInModalSpace(const int &vertexID, const Eigen::Vector3d &impulse, Eigen::VectorXd &forceInModalSpace);
        // Get u from q by doing u = Uq, where U is the modal matrix.
        void GetVolumeVertexDisplacement(const Eigen::VectorXd &q, Eigen::VectorXd &u); 
        void GetVolumeVertexModeValues(const int &modeIndex, Eigen::VectorXd &modeValues); 
        void SetODESolverTime(const REAL &time);
        void Initialize(const REAL &ODEStepSize, const std::string &modeFile, std::shared_ptr<ModalMaterial> materialPtr);
        void InitializeModalODESolvers(std::shared_ptr<ModalMaterial> materialPtr);
        void InitializeSparseModalEncoder(); 
        // remove all non-surface mode shapes and permute the rows of modal
        // matrix
        void CullNonSurfaceModeShapes(std::shared_ptr<TetMeshIndexToSurfaceMesh> idMapPtr, std::shared_ptr<TriangleMesh<REAL> > meshPtr); 
        void CullHighFrequencyModes(const int &modesToKeep); 
        void CullHighFrequencyModes(const REAL &frequenciesToKeep); 

    friend std::ostream &operator <<(std::ostream &os, const ModalAnalysisObject &object); 
}; 

#endif
