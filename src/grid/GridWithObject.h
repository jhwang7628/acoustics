#ifndef GRIDWITHOBJECT
#define GRIDWITHOBJECT

#include <string>
#include <iostream> 

#include "Grid.h"
#include "parser/Parser.h"
#include <distancefield/closestPointField.h> 
#include <distancefield/FieldBuilder.h>

#include <linearalgebra/Vector3.hpp> 

// overwrites gradient computation in uniform grid because the presence of the
// object needs to be handled. 
class UniformGridWithObject : public UniformGrid
{

    public: 
        enum CellType : bool { FLUID=false, SOLID=true };
        //typedef typename UniformGrid<T>::CELL_GRADIENT_COMPONENT CELL_GRADIENT_COMPONENT; 
        //typedef UniformGrid<T>::CELL_GRADIENT_COMPONENT CELL_GRADIENT_COMPONENT; 

    protected: 
        std::string                             meshName_; 
        std::shared_ptr<ClosestPointField>      distanceField_; 
        std::shared_ptr<TriangleMesh<double>>   mesh_; 
        std::shared_ptr<Parser>                 parser_; 
        Parser::ImpulseResponseParms            solverParameters_; 

        std::vector<CellType>                   isBoundary_; 

        double                                  distanceTolerance_; 
        bool                                    initialized_; 

    public:
        UniformGridWithObject(const Eigen::Vector3d &minBound, const Eigen::Vector3d &maxBound, const Eigen::Vector3i &cellCount) : 
            UniformGrid(minBound, maxBound, cellCount), 
            initialized_(false) { }

        UniformGridWithObject() { } 

        void Reinitialize(const std::string &configFile); 
        void ClassifyCells();

        // same as gradient in UniformGrid class except with object handling.
        // For function usage please refer to that class.
        virtual void CellCenteredDataGradient( const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &gradientData , const CELL_GRADIENT_COMPONENT &component=ALL); 

};


#endif
