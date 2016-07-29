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
        // i.e., rows: #frames; 1st index
        //       cols: #objects;  2nd index
        typedef std::vector<REAL>               VecReal; 
        typedef std::vector<int>                VecInt; 
        typedef std::vector<Point3<REAL> >      VecPoint3; 
        typedef std::vector<Vector3<REAL> >     VecVector3; 
        typedef std::vector<Quaternion<REAL> >  VecQuaternion; 
        typedef std::vector<VecPoint3>          MatPoint3; 
        typedef std::vector<VecVector3>         MatVector3; 
        typedef std::vector<VecQuaternion>      MatQuaternion; 

    private: 
        int             _totalFrames; 
        int             _totalObjects; 

        VecReal         _timesteps; 
        MatPoint3       _positions; 
        MatQuaternion   _rotations; 
        MatVector3      _velocity; 
        MatVector3      _angularVelocity; 
        MatVector3      _acceleration; 
        MatVector3      _angularAcceleration; 

    public: 
        RigidObjDispReader()
            : _totalFrames(0), _totalObjects(0)
        { }

        inline int N_Frames(){return _totalFrames;}
        inline int N_Objects(){return _totalObjects;}

        // getters
        inline VecReal &Timesteps(){return _timesteps;}
        inline const VecReal &Timesteps() const {return _timesteps;}
        inline MatPoint3 &Positions(){return _positions;}
        inline const MatPoint3 &Positions() const {return _positions;}
        inline MatQuaternion &Rotations(){return _rotations;}
        inline const MatQuaternion &Rotations() const {return _rotations;}

        void ReadDisplacement(const std::string &filename); 
        void ReadAllKinematics(const std::string &fileDisplacement, const std::string &fileVelocity, const std::string &fileAcceleration);  
}; 

#endif
