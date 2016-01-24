#ifndef GRIDWITHOBJECT
#define GRIDWITHOBJECT

#include <string>
#include <iostream> 

#include "Grid.h"
#include "parser/Parser.h"
#include <distancefield/closestPointField.h> 
#include <distancefield/FieldBuilder.h>



// overwrites gradient computation in uniform grid because the presence of the
// object needs to be handled. 
template <class T> 
class UniformGridWithObject : public UniformGrid<T>
{

    protected: 
        std::string                             meshName_; 
        std::shared_ptr<ClosestPointField>      distanceField_; 
        std::shared_ptr<TriangleMesh<double>>   mesh_; 
        std::shared_ptr<Parser>                 parser_; 

    public:
        UniformGridWithObject(const Eigen::Vector3d &minBound, const Eigen::Vector3d &maxBound, const Eigen::Vector3i &cellCount) : 
            UniformGrid<T>(minBound, maxBound, cellCount) { }
        UniformGridWithObject() { } 

        void Reinitialize(const std::string configFile)
        {
            std::cout << "config file : " << configFile << std::endl;

            Parser::AcousticTransferParms parms; 

            parser_.reset(Parser::buildParser( configFile )); 
            if ( !parser_ ) throw std::runtime_error("**ERROR** Could not build parser from "+configFile);

            mesh_.reset(parser_->getMesh()); 
            if ( !mesh_ ) throw std::runtime_error("**ERROR** Could not build mesh");

            parms = parser_->getAcousticTransferParms();

            distanceField_.reset(DistanceFieldBuilder::BuildSignedClosestPointField(parser_->getMeshFileName().c_str(), parms._sdfResolution, parms._sdfFilePrefix.c_str())); 

            if (!distanceField_) throw std::runtime_error("**ERROR** Could not construct distance field"); 


        }



};


#endif
