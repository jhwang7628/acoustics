#include <iostream> 
#include <fstream> 
#include <io/RigidObjDispReader.h> 

//##############################################################################
//##############################################################################
void RigidObjDispReader::
ReadDisplacement(const std::string &filename) 
{
    std::cout << "Read displacement from rigidsim data : " << filename << "\n"; 
    std::ifstream inFile(filename, std::ios::in | std::ios::binary); 
    if (!inFile) throw std::runtime_error("**ERROR** Can't open file: " + filename); 

    REAL time; 
    if (!inFile.eof())
        inFile.read((char*)&time, sizeof(REAL));
    else 
        throw std::runtime_error("**WARNING** No time step found in the file"); 

    _totalObjects = 0; 
    _totalFrames = 0;
    while(true)
    {
        int objID = -1; 
        int totalObjects = 0; 
        VecPoint3 framePositions; 
        VecQuaternion frameRotations; 
        inFile.read((char*)&objID, sizeof(int)); 
        while(objID != -1) // end flag for the time step
        {
            Point3<REAL>        displacement; 
            Quaternion<REAL>    rotation; 
            inFile.read((char*)&displacement, sizeof(Point3<REAL>)); 
            inFile.read((char*)&rotation, sizeof(Quaternion<REAL>)); 
            framePositions.push_back(displacement); 
            frameRotations.push_back(rotation); 
            totalObjects ++;
            inFile.read((char*)&objID, sizeof(int)); 
        }
        _timesteps.push_back(time); 
        _positions.push_back(framePositions); 
        _rotations.push_back(frameRotations); 
        _totalFrames ++; 
        if (totalObjects != 0) _totalObjects = totalObjects; 

        if (!inFile.eof())
            inFile.read((char*)&time, sizeof(REAL));
        else 
            break; 
    }
    std::cout << " Frames read : " << _totalFrames << std::endl;
    std::cout << " Objects read: " << _totalObjects << std::endl;
}

//##############################################################################
//##############################################################################
void RigidObjDispReader::
ReadAllKinematics(const std::string &fileDisplacement, const std::string &fileVelocity, const std::string &fileAcceleration) 
{
    std::cout << "read all kinematics\n"; 
}
