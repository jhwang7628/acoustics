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
        typedef unsigned char CellType; 
        // bit flags : http://www.cplusplus.com/forum/general/1590/
        // [1|2|3|4|5|6|7|8]: 
        //   1: whether its solid cell: enclosed by object using sdf
        //   2: whether its interface cell: solid cell that has at least one
        //      neighbour that is not solid 
        //   3: if cell has solid on left  in x-direction, 
        //   4: if cell has solid on right in x-direction.
        //   5: if cell has solid on left  in y-direction, 
        //   6: if cell has solid on right in y-direction.
        //   7: if cell has solid on left  in z-direction, 
        //   8: if cell has solid on right in z-direction.
        enum InterfacialInfo {  IS_SOLID     =0x80,
                                IS_INTERFACE =0x40,
                                X_SOLID_ON_LEFT =0x20, 
                                X_SOLID_ON_RIGHT=0x10, 
                                Y_SOLID_ON_LEFT =0x08, 
                                Y_SOLID_ON_RIGHT=0x04, 
                                Z_SOLID_ON_LEFT =0x02, 
                                Z_SOLID_ON_RIGHT=0x01, }; 

    protected: 
        std::string                             meshName_; 
        std::shared_ptr<ClosestPointField>      distanceField_; 
        std::shared_ptr<TriangleMesh<double>>   mesh_; 
        std::shared_ptr<Parser>                 parser_; 
        Parser::ImpulseResponseParms            solverParameters_; 

        std::vector<CellType>                   cellTypes_;  
        double                                  distanceTolerance_; 
        bool                                    initialized_; 

    public:
        UniformGridWithObject(const Eigen::Vector3d &minBound, const Eigen::Vector3d &maxBound, const Eigen::Vector3i &cellCount) : 
            UniformGrid(minBound, maxBound, cellCount), 
            initialized_(false) { }

        UniformGridWithObject() { } 

        void Reinitialize(const std::string &configFile); 
        void ClassifyCells();

        void WriteCellTypes(const std::string &filename); 

        void InterfacialGradientSmoothing();

        // same as gradient in UniformGrid class except with object handling.
        // For function usage please refer to that class.
        virtual void CellCenteredDataGradient( const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &gradientData , const CELL_GRADIENT_COMPONENT &component=ALL); 

};


#endif
