#ifndef MODAL_ANALYSIS_H 
#define MODAL_ANALYSIS_H 
#include <memory>
#include <sndgen/RigidModal.h> 

//##############################################################################
// Class that handles the modal objects 
//##############################################################################
class ModalAnalysis
{
    public: 
        typedef std::shared_ptr<RigidModal> RigidModalPtr; 

    private: 
        std::string                 _configFile; 
        std::vector<RigidModalPtr>  _rigidModals; 

    public: 
        ModalAnalysis(const std::string &configFile) 
            : _configFile(configFile)
        {}

        inline RigidModalPtr GetModal(const int &index){return _rigidModals.at(index);}
        inline const RigidModalPtr &GetModal(const int &index) const {return _rigidModals.at(index);}
        void BuildModalModelsFromFile(); 
}; 

#endif
