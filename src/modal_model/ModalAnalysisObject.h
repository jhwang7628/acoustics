#ifndef MODAL_ANALYSIS_OBJECT_H 
#define MODAL_ANALYSIS_OBJECT_H 
#include <iostream>
#include <memory>
#include <deformable/ModeData.h> 
#include <Eigen/Dense>

//##############################################################################
// Class that handles the modal object
//##############################################################################
class ModalAnalysisObject
{
    private: 
        bool        _modeFileSet; 
        std::string _modeFile; 

        Eigen::VectorXd _eigenValues; // eigenvalues produced by modal analysis (omega squared)
        Eigen::MatrixXd _eigenVectors; 

    public: 
        ModalAnalysisObject()
            : _modeFileSet(false)
        {} 
        ModalAnalysisObject(const std::string &modeFile) 
            : _modeFile(modeFile), _modeFileSet(true)
        {}

        inline int N_Modes() const {return _eigenVectors.cols();}
        inline int N_DOF() const {return _eigenVectors.rows();}
        inline int N_vertices() const {return N_DOF()/3;}
        inline void SetModeFile(const std::string &modeFile){_modeFile = modeFile; _modeFileSet = true;}
        void ReadModeFromFile(); 
        // Get U^T f, where U is the eigenvector (modal) matrix. Input vertexID should be in zero-based. 
        void GetForceInModalSpace(const int &vertexID, const Eigen::Vector3d &impulse, Eigen::VectorXd &forceInModalSpace);
        // Get u from q by doing u = Uq, where U is the modal matrix.
        void GetVertexDisplacement(const Eigen::VectorXd &q, Eigen::VectorXd &u); 

    friend std::ostream &operator <<(std::ostream &os, const ModalAnalysisObject &object); 
}; 

#endif
