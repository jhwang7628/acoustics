#ifndef PRESSURE_SOURCE_H 
#define PRESSURE_SOURCE_H 

#include <TYPES.h> 
#include <config.h>
#include <Eigen/Dense> 

#include <config.h> 
#include <linearalgebra/Vector3.hpp> 
#include <wavesolver/Source.h> 

//##############################################################################
// Injected pressure source
//##############################################################################
class PressureSource : public Source 
{
    public: 
        PressureSource(){} 

        // normal would not be used only exists for consistency with
        // vibrational source. can pass in, e.g., position twice
        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle=-1)=0; 
        virtual REAL Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time)=0; 
        virtual Vector3d Evaluate(const int &vertexID, const REAL &time)=0;
        virtual void PrintSourceInfo(std::ostream &os)const=0; 
};

#endif
