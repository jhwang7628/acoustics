#ifndef MODAL_ANALYSIS_OBJECT_H 
#define MODAL_ANALYSIS_OBJECT_H 
#include <iostream>
#include <memory>
#include <deformable/ModeData.h> 
#include <modal_model/ModalODESolver.h> 
#include <modal_model/ModalMaterial.h>
#include <Eigen/Dense>

//##############################################################################
// Class that handles the modal object
//
// Note: 
//  eigenvectors are based on tet meshes
//##############################################################################
class ModalAnalysisObject
{
    public: 
        typedef std::shared_ptr<ModalODESolver> ModalODESolverPtr; 

    private: 
        bool        _modeFileSet; 
        std::string _modeFile; 

        Eigen::VectorXd _eigenValues; // eigenvalues produced by modal analysis (omega squared)
        Eigen::MatrixXd _eigenVectors; 

        std::vector<ModalODESolverPtr> _modalODESolvers; 
        REAL        _timeStepSize; 

    public: 
        ModalAnalysisObject()
            : _modeFileSet(false)
        {} 
        ModalAnalysisObject(const std::string &modeFile) 
            : _modeFileSet(true), _modeFile(modeFile)
        {}

        inline int N_Modes() const {return _eigenVectors.cols();}
        inline int N_DOF() const {return _eigenVectors.rows();}
        inline int N_vertices() const {return N_DOF()/3;}
        inline void SetModeFile(const std::string &modeFile){_modeFile = modeFile; _modeFileSet = true;}
        void ReadModeFromFile(); 
        // Get U^T f, where U is the eigenvector (modal) matrix. Input vertexID should be in zero-based, tet mesh id. 
        void GetForceInModalSpace(const int &vertexID, const Eigen::Vector3d &impulse, Eigen::VectorXd &forceInModalSpace);
        // Get u from q by doing u = Uq, where U is the modal matrix.
        void GetVolumeVertexDisplacement(const Eigen::VectorXd &q, Eigen::VectorXd &u); 
        void GetVolumeVertexModeValues(const int &modeIndex, Eigen::VectorXd &modeValues); 
        void Initialize(const std::string &modeFile, std::shared_ptr<ModalMaterial> materialPtr);
        void InitializeModalODESolvers(std::shared_ptr<ModalMaterial> materialPtr);
        void StepAllModalODESolvers(const int &vertexID, const Eigen::Vector3d &impulse); 

    friend std::ostream &operator <<(std::ostream &os, const ModalAnalysisObject &object); 
}; 

#endif
