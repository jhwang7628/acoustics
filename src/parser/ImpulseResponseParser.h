#ifndef IMPULSE_RESPONSE_PARSER_H 
#define IMPULSE_RESPONSE_PARSER_H 

#include <string> 
#include <vector> 
#include <TYPES.h> 
#include <tinyxml/tinyxml.h> 
#include <tinyxml/tinystr.h> 
#include <parser/Parser.h>
#include <wavesolver/VolumetricSource.h> 
#include <wavesolver/FDTD_Objects.h> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h> 

//##############################################################################
// Parameters for impulse response wave solve
//##############################################################################
class ImpulseResponseParser : public Parser
{
    public: 
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

            std::vector<VolumetricSource>   _sources; 

            bool                            _useMesh; 
            bool                            _cornellBoxBoundaryCondition; // if on then PML is only on one face. all the other boundaries for the box are hard-wall
            bool                            _useGhostCellBoundary; // if on then ghost cell will be used, otherwise rasterized boundary

            REAL                            _f; 

            void SetSources(const std::vector<VolumetricSource> &sources) { _sources = sources; } 
        };


    private: 


    public: 
        ImpulseResponseParser(){}
        ImpulseResponseParser(const std::string &documentName)
        {
            TiXmlDocument *document = new TiXmlDocument(); 
            document->LoadFile(documentName.c_str()); 
            SetDocument(document); 
        }

        static std::vector<VolumetricSource> QueryVolumetricSource(TiXmlNode *document, Parser *parser, const std::string &path, const REAL &soundSpeed); 
        inline void SetDocument(TiXmlDocument *document){this->document = document;}
        inline TiXmlDocument *GetDocument(){return document;}
        // Fetch parameters for impulse response computation.
        ImpulseResponseParms Parse(); 
        void GetObjects(const std::string &inputElement, FDTD_Objects &objects); 
        //void getMeshes( const std::string &inputElement, vector<TriangleMesh<REAL> *> &meshes, 
        //                vector<string> &rigidFilePrefixes, vector<string> &sdfFileNames, 
        //                IntArray &sdfResolutions, vector<string> &pulseModelFileNames, 
        //                FloatArray &densities, map<std::string, int> &meshIDMap, bool compressedFields ); 



}; 

#endif 
