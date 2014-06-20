// Parser.h: Interface for the Parser class
//
//////////////////////////////////////////////////////////////////////

#ifndef PARSER_H
#define PARSER_H

#include <tinyxml/tinyxml.h>
#include <tinyxml/tinystr.h>

#include <TYPES.h>

#include <geometry/ClosestPointMesh.h>
#include <geometry/GTS_TriMesh.h>
#include <geometry/RigidMesh.h>
#include <geometry/TriangleMesh.hpp>

#include <transfer/RadialApproximation.h>

#include <utils/IO.h>
#include <utils/STLUtil.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

// A helper macro
#define FOR_CHILD_ELEMENTS( childVar, parentVar, childValue ) \
    for( TiXmlElement* childVar = parentVar->FirstChildElement(childValue); \
            childVar; \
            childVar = childVar->NextSiblingElement(childValue) )

//////////////////////////////////////////////////////////////////////
// Parser class
//
// A class which can parse input files and produce triangle meshes,
// basis information, etc. for a given example
//////////////////////////////////////////////////////////////////////
class Parser {
    public:
        static Parser *buildParser(std::string filename);

        // Construct a triangle mesh from information stored
        // in this document
        TriangleMesh<REAL> *getMesh();

        // Mesh file name
        std::string getMeshFileName();

#if 0
        struct EquivalentSourceParms {
            // The acceleration pulse function
            Function               _acceleration;

            int                    _multipoleTerms;

            // Multipole pulse time scale
            REAL                   _timeScale;

            // Width between pulses
            REAL                   _pulseInterval;

            // Integration time step
            REAL                   _timeStep;

            // Speed of sound and density
            REAL                   _c;
            REAL                   _density;
        };
#endif

#if 0
        struct SourceTrainingParms {
            // The acceleration pulse function
            Function               _acceleration;

            // Multipole pulse time scale
            REAL                   _timeScale;

            // Integration time step
            REAL                   _timeStep;

            // Speed of sound and density
            REAL                   _c;
            REAL                   _density;

            /////////////////////////////////////////////////////
            // SDF parms
            /////////////////////////////////////////////////////
            int                    _sdfResolution;

            std::string            _sdfFilePrefix;

            int                    _pointsPerVoxel;

            /////////////////////////////////////////////////////
            // Training parms
            /////////////////////////////////////////////////////
            int                    _sourcesPerIteration;

            int                    _candidatesPerIteration;

            int                    _maxSourcesPerBatch;

            REAL                   _errorTolerance;

            void printParms( ostream &os ) const
            {
                os << "Training parameters:" << endl;
                os << SDUMP( _timeScale ) << endl;
                os << SDUMP( _timeStep ) << endl;
                os << SDUMP( _c ) << endl;
                os << SDUMP( _density ) << endl;
                os << SDUMP( _sdfResolution ) << endl;
                os << SDUMP( _sdfFilePrefix ) << endl;
                os << SDUMP( _pointsPerVoxel ) << endl;
                os << SDUMP( _sourcesPerIteration ) << endl;
                os << SDUMP( _candidatesPerIteration ) << endl;
                os << SDUMP( _maxSourcesPerBatch ) << endl;
                os << SDUMP( _errorTolerance ) << endl;
                os << endl;
            }

        };
#endif

#if 0
        struct SphericalImpactParms {
            /////////////////////////////////////////////////////
            // Speed of sound and density
            /////////////////////////////////////////////////////
            REAL                   _c;
            REAL                   _density;

            /////////////////////////////////////////////////////
            // SDF parms
            /////////////////////////////////////////////////////
            int                    _sdfResolution;

            std::string            _sdfFilePrefix;

            /////////////////////////////////////////////////////
            // Domain resolution and size
            /////////////////////////////////////////////////////
            int                    _gridResolution;

            // How much to scale the bounding box of the object by
            REAL                   _gridScale;

            int                    _timeStepFrequency;

            /////////////////////////////////////////////////////
            // Boundary information
            /////////////////////////////////////////////////////

            // Radius of "sphere"
            REAL                   _sphereRadius;

            // Impact velocity
            REAL                   _impactSpeed;

            // Material parameters (Young's modulus and Poisson ratio)
            REAL                   _boundary_E;
            REAL                   _boundary_nu;

            // Sphere mass
            REAL                   _boundaryDensity;

            /////////////////////////////////////////////////////
            // Sound output parameters
            /////////////////////////////////////////////////////
            VEC3F                  _listeningPosition;
            Vector3Array           _listeningPositions;

            std::string            _outputFile;

        };
#endif

        /////////////////////////////////////////////////////
        // Parameters for acceleration pulse precomputation
        /////////////////////////////////////////////////////
        struct AcousticTransferParms {
            /////////////////////////////////////////////////////
            // Speed of sound and density
            /////////////////////////////////////////////////////
            REAL                   _c;
            REAL                   _density;

            /////////////////////////////////////////////////////
            // SDF parms
            /////////////////////////////////////////////////////
            int                    _sdfResolution;

            std::string            _sdfFilePrefix;

            /////////////////////////////////////////////////////
            // Domain resolution and size
            /////////////////////////////////////////////////////
            int                    _gridResolution;

            // How much to scale the bounding box of the object by
            REAL                   _gridScale;

            // Whether or not to use a fixed grid cell size
            bool                   _fixedCellSize;
            REAL                   _cellSize;

            int                    _timeStepFrequency;
            int                    _subSteps;

            /////////////////////////////////////////////////////
            // For computing a multi-term expansion
            /////////////////////////////////////////////////////
            int                    _nbar;
            REAL                   _radiusMultipole;
            std::string            _modeDataFile;

            /////////////////////////////////////////////////////
            // Rigid mesh information
            /////////////////////////////////////////////////////
            std::string            _rigidPrefix;
            REAL                   _rigidDensity;

            /////////////////////////////////////////////////////
            // Output parameters
            /////////////////////////////////////////////////////
            std::string            _outputFile;

            std::string            _multipoleOutputFile;

        };

        /////////////////////////////////////////////////////
        // Parameters for acceleration pulse precomputation
        /////////////////////////////////////////////////////
        struct AccelerationPulseParms {
            /////////////////////////////////////////////////////
            // Speed of sound and density
            /////////////////////////////////////////////////////
            REAL                   _c;
            REAL                   _density;

            /////////////////////////////////////////////////////
            // Time scale of the pulse to precompute
            /////////////////////////////////////////////////////
            REAL                   _pulseTimeScale;

            /////////////////////////////////////////////////////
            // SDF parms
            /////////////////////////////////////////////////////
            int                    _sdfResolution;

            std::string            _sdfFilePrefix;

            /////////////////////////////////////////////////////
            // Domain resolution and size
            /////////////////////////////////////////////////////
            int                    _gridResolution;

            // How much to scale the bounding box of the object by
            REAL                   _gridScale;

            // Whether or not to use a fixed grid cell size
            bool                   _fixedCellSize;
            REAL                   _cellSize;

            int                    _timeStepFrequency;
            int                    _subSteps;

            /////////////////////////////////////////////////////
            // Radial approximation parameters
            /////////////////////////////////////////////////////
            int                    _listeningResolution;
            REAL                   _listeningRadius;

            /////////////////////////////////////////////////////
            // For computing a multi-term expansion
            /////////////////////////////////////////////////////
            int                    _numShells;
            int                    _numTerms;
            REAL                   _radiusMultiplier;

            /////////////////////////////////////////////////////
            // Rigid mesh information
            /////////////////////////////////////////////////////
            std::string            _rigidPrefix;
            REAL                   _rigidDensity;

            /////////////////////////////////////////////////////
            // If we want to do wavelet compression
            /////////////////////////////////////////////////////
            REAL                   _compressionTolerance;

            std::string            _compressedOutputFile;

            bool                   _symmetrizeField;
            int                    _symmetricResolution;

            /////////////////////////////////////////////////////
            // Use this if the field is only represented in
            // a specific set of directions
            /////////////////////////////////////////////////////
            bool                   _useDirectionSet;
            std::string            _directionSetFile;

            /////////////////////////////////////////////////////
            // Output parameters
            /////////////////////////////////////////////////////
            std::string            _outputFile;

            std::string            _multiTermOutputFile;

        };

#if 0
        struct ScatterObject {
            std::string          _fileName;

            REAL                 _scale;
            VEC3F                _translation;

            VEC3F                _accelerationDirection;
            REAL                 _accelerationMagnitude;

            REAL                 _absorptionCoefficient;

            /////////////////////////////////////////////////////
            // Required attributes for field computation if this
            // object has a non-zero acceleration magnitude
            /////////////////////////////////////////////////////

            /////////////////////////////////////////////////////
            // SDF parms
            /////////////////////////////////////////////////////
            int                    _sdfResolution;

            std::string            _sdfFilePrefix;

            /////////////////////////////////////////////////////
            // Domain resolution and size
            /////////////////////////////////////////////////////
            int                    _gridResolution;

            // How much to scale the bounding box of the object by
            REAL                   _gridScale;

            int                    _timeStepFrequency;

            // Substeps to actually take (but not store)
            int                    _subSteps;

            /////////////////////////////////////////////////////
            // Boundary information
            /////////////////////////////////////////////////////

            // Radius of "sphere"
            REAL                   _sphereRadius;

            // Impact velocity
            REAL                   _impactSpeed;

            // Material parameters (Young's modulus and Poisson ratio)
            REAL                   _boundary_E;
            REAL                   _boundary_nu;

            // Sphere mass
            REAL                   _boundaryDensity;
        };
#endif

#if 0
        struct ScatteringTestParms {
            /////////////////////////////////////////////////////
            // Speed of sound and density
            /////////////////////////////////////////////////////
            REAL                   _c;
            REAL                   _density;

            /////////////////////////////////////////////////////
            // Read files with input frequencies and their
            // coefficients
            /////////////////////////////////////////////////////
            std::string            _frequencyFile;

            /////////////////////////////////////////////////////
            // Object file names, acceleration magnitudes, and
            // directions for objects in the scene
            /////////////////////////////////////////////////////
            std::vector<ScatterObject> _scatterObjects;

            /////////////////////////////////////////////////////
            // Output file and listening positions
            /////////////////////////////////////////////////////
            std::string            _outputPrefix;
            std::string            _spectrumOutputPrefix;

            Vector3Array           _listeningPositions;

            /////////////////////////////////////////////////////
            // Output data files
            /////////////////////////////////////////////////////
            std::string            _responseFile;

            /////////////////////////////////////////////////////
            // Finite difference grid geometry
            /////////////////////////////////////////////////////
            VEC3F                  _FD_min_bound;
            VEC3F                  _FD_max_bound;
            int                    _FD_grid_divisions;

            /////////////////////////////////////////////////////
            // Finite difference time stepping
            /////////////////////////////////////////////////////
            int                    _FD_timeStepFrequency;

            // Substeps to actually take (but not store)
            int                    _FD_subSteps;

            std::string            _FD_outputFile;

            // Impact time scale for FD simulation
            REAL                   _FD_impactTimeScale;

        };
#endif

        struct SceneObjectParms {
            SceneObjectParms()
            {
            }

            void clear()
            {
                clearVectorContents( _distanceMeshes );
                clearVectorContents( _curvatureMeshes );
                clearVectorContents( _rigidMeshes );
                clearVectorContents( _meshes );
                clearVectorContents( _pulseModels );
            }

            // All triangle meshes in scene
            std::vector<TriangleMesh<REAL> *>  _meshes;

            // Structure built around each mesh for distance queries
            std::vector<ClosestPointMesh *>    _distanceMeshes;

            // Structure built around each mesh for curvature queries
            std::vector<GTS_TriMesh *>         _curvatureMeshes;

            // Structure built around each mesh to store rigid body information
            std::vector<RigidMesh *>           _rigidMeshes;

            // Acceleration pulse approximations for each mesh
            std::vector<RadialApproximation::AccelerationSet *>
                _pulseModels;

            std::vector<CompactRadialApproximation *>
                _compactPulseModels;

            // Each object in the scene has a mesh ID which refers to the
            // meshes in the set above
            IntArray                           _objectMeshIDs;

            // If we are using Wavelet-compressed PAN, enable this parameter
            // to truncate high-frequency wavelet coefficients and match
            // PAN sampling rates (approximately) to the output sampling rate
            bool                               _matchWaveletSamplingRates;

            // Proxy parameters
            bool                               _useProxies;

            ObjectMeasure                      _proxyMeasure;

            std::string                        _proxyPathPrefix;
            std::string                        _proxyFilePrefix;

            REAL                               _proxyMinScale;
            int                                _proxyIncrements;

            // Symmetrized proxies
            bool                               _useSymmetricProxies;
            std::string                        _symmetricProxyPrefix;

            // In case we want to use a proxy model in which PAN fields
            // are discretized over a direction set (rather than over
            // a uniform angular distribution)
            bool                               _useDirectionSetProxies;
            std::string                        _proxyDirectionSetFile;
            std::string                        _proxyDirectionSetSuffix;

            // Parameters for impulse time randomization
            bool                               _randomizeImpulseTimes;
            REAL                               _simTimeStep;

            // Optional time scaling argument
            REAL                               _collisionTimeScale;

        };

#if 0
        // Fetches equivalent source parameters from the XML document
        EquivalentSourceParms getEquivalentSourceParms();
#endif

#if 0
        // Fetches parameters for source position training
        SourceTrainingParms getSourceTrainingParms();
#endif

#if 0
        // Fetches parameters for faking a sphere-like collision (with an
        // arbitrary mesh)
        SphericalImpactParms getSphereImpactParms();
#endif

        // Fetches parameters for acceleration pulse precomputation
        AcousticTransferParms getAcousticTransferParms();
        AccelerationPulseParms getAccelerationPulseParms();

#if 0
        // Fetches parameters for setting up the input to an acoustic
        // scattering calculation due to rigid body acceleration
        ScatteringTestParms getScatteringTestParms();
#endif

        // Fetches scene objects
        SceneObjectParms getSceneObjectParms( bool loadPANData,
                                              bool compressedFieldse,
                                              bool symmetricFieldse );

        // If the attribute DNE, it will print out an error message, prompt
        // for input, and return "ERROR" for the value.
        static std::string queryRequiredAttr( TiXmlElement* node, std::string attr );

        // If the attribute DNE, it will return the given defaultValue
        static std::string queryOptionalAttr( TiXmlElement* node, std::string attr,
                std::string defaultValue );

        // Lame implementation of simple xpath. Uses the first child if there are
        // multiple elements matching a given name.
        TiXmlNode* getNodeByPath( std::string pathStr );

        // Convenience function. Example:
        // return queryRequiredAttr( string("shell/farfieldsound"),
        // string("outputDir") );
        std::string queryRequiredAttr( std::string path, std::string attr );
        std::string queryOptionalAttr( std::string path, std::string attr,
                std::string defaultVal );

        std::string queryOptionalAttr( const char* path, const char* attr,
                const char* defaultVal )
        {
            return queryOptionalAttr( std::string(path), std::string(attr),
                    std::string(defaultVal) );
        }

        std::string queryRequiredAttr( const char* path, const char* attr )
        {
            return queryRequiredAttr( std::string(path), std::string(attr) );
        }

        REAL queryOptionalReal( const char* path, const char* attr,
                const char* defaultVal )
        {
            return atof( queryOptionalAttr( path, attr, defaultVal ).c_str() );
        }
        REAL queryRequiredReal( const char* path, const char* attr )
        {
            return atof( queryRequiredAttr( path, attr ).c_str() );
        }

        REAL queryOptionalInt( const char* path, const char* attr,
                const char* defaultVal )
        {
            return atoi( queryOptionalAttr( path, attr, defaultVal ).c_str() );
        }
        REAL queryRequiredInt( const char* path, const char* attr )
        {
            return atoi( queryRequiredAttr( path, attr ).c_str() );
        }

        // If you want to walk the DOM on your own, use this
        TiXmlDocument* getDocument() { return document; }

        virtual ~Parser();

    private:
        // Gets a list of listening positions from the given input element
        void getListeningPositions( const char *inputElement,
                                    Vector3Array &positions );

        // Reads a mesh list
        //
        // Helper function for getSceneObjectParms
        void getMeshes( const char *inputElement,
                        std::vector<TriangleMesh<REAL> *> &meshes,
                        std::vector<string> &rigidFilePrefixes,
                        std::vector<string> &sdfFileNames,
                        IntArray &sdfResolutions,
                        std::vector<string> &pulseModelFileNames,
                        FloatArray &densities,
                        std::map<std::string, int> &meshIDMap,
                        bool compressedFields );

        // Gets the mesh ID for each object in a scene
        void getObjectMeshIDs( const char *inputElement,
                               std::map<std::string, int> &meshIDMap,
                               IntArray &objectMeshIDs );

#if 0
        // Gets a list of scatter objects from the given input element
        void getScatterObjects( const char *inputElement,
                std::vector<ScatterObject> &objects );
#endif

        // Construct a new parser given a tinyxml document.  We do
        // this privately so that we can make sure a file exists
        // before constructing the parser
        Parser(TiXmlDocument *document);

        TiXmlDocument *document;

        // Kept mainly for helpful error reporting
        std::string filename;
};

#endif
