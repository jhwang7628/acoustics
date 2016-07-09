#ifndef RIGID_OBJ_DISP_READER_H 
#define RIGID_OBJ_DISP_READER_H 
#include <vector> 
#include <config.h> 
#include <geometry/Point3.hpp> 
#include <linearalgebra/Vector3.hpp> 
#include <linearalgebra/Quaternion.hpp> 

//##############################################################################
// Class for reading files saved by class RigidObjDispRecorder.
// Want Structures of Arrays so interfacing with other class is easier. 
//##############################################################################
class RigidObjDispReader
{
    public: 
        // all of the matrix data is of dimension #objects * #frames
        // i.e., rows: #objects; 1st index
        //       cols: #frames;  2nd index
        typedef std::vector<std::vector<Point3<REAL> > >        MatPoint3; 
        typedef std::vector<std::vector<Vector3<REAL> > >       MatVector3; 
        typedef std::vector<std::vector<Quaternion<REAL> > >    MatQuaternion; 

    private: 
        int             _totalFrames; 
        int             _totalObjects; 

        MatPoint3       _positions; 
        MatQuaternion   _rotations; 
        MatVector3      _velocity; 
        MatVector3      _angularVelocity; 
        MatVector3      _acceleration; 
        MatVector3      _angularAcceleration; 

    public: 
        void ReadDisplacement(const std::string &filename); 
        void ReadAllKinematics(const std::string &fileDisplacement, const std::string &fileVelocity, const std::string &fileAcceleration);  
}; 

#endif
