#ifndef IMPULSE_RESPONSE_PARSER_H 
#define IMPULSE_RESPONSE_PARSER_H 

#include <string> 
#include <vector> 
#include <TYPES.h> 
#include <tinyxml/tinyxml.h> 
#include <tinyxml/tinystr.h> 
#include <parser/Parser.h>
#include <modal_model/ModalMaterial.h>
#include <modal_model/ModalMaterialList.h>
#include <wavesolver/PressureSource.h> 
#include <wavesolver/FDTD_Objects.h> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h> 
#include <wavesolver/PML_WaveSolver_Settings.h> 

//##############################################################################
// Parameters for impulse response wave solve
//##############################################################################
class ImpulseResponseParser : public Parser
{
    public: 
        // this struct is kept for backward compatibility,
        // tools/impulse-response will need it I think.
        struct ImpulseResponseParms 
        {
            // Speed of sound and density
            REAL                            _c;
            REAL                            _density;

            // SDF parms
            int                             _sdfResolution;
            std::string                     _sdfFilePrefix;

            // Domain resolution and         size
            int                             _gridResolution;

            // How much to scale the         bounding box of the object by
            REAL                            _gridScale;

            // Fix cell size (optional parameter)
            REAL                            _cellSize;

            int                             _timeStepFrequency;
            int                             _subSteps;

            REAL                            _stopTime; 

            // Output parameters
            std::string                     _outputPattern; 

            std::string                     _listeningFile; 

            //std::vector<VolumetricSource>   _sources; 

            bool                            _useMesh; 
            bool                            _cornellBoxBoundaryCondition; // if on then PML is only on one face. all the other boundaries for the box are hard-wall
            bool                            _useGhostCellBoundary; // if on then ghost cell will be used, otherwise rasterized boundary

            REAL                            _f; 

            //void SetSources(const std::vector<VolumetricSource> &sources) { _sources = sources; } 
        };

    private:
        TiXmlDocument _document; 

    public: 
        ImpulseResponseParser(){}
        ImpulseResponseParser(const std::string &documentName)
        {
            _document.LoadFile(documentName.c_str()); 
        }

        inline void SetDocument(const TiXmlDocument &document){this->_document = document;}
        inline TiXmlDocument &GetDocument(){return _document;}
        void GetObjects(std::shared_ptr<FDTD_Objects> &objects); 
        //void GetSources(std::shared_ptr<PML_WaveSolver_Settings> settings);
        void GetSolverSettings(std::shared_ptr<PML_WaveSolver_Settings> &settings); 
        void GetPressureSources(const REAL &soundSpeed, std::vector<PressureSourcePtr> &pressureSources); 
        void GetListeningPoints(Vector3Array &listeningPoints); 
        void GetModalMaterials(ModalMaterialList &modalMaterials); 
}; 

#endif 
