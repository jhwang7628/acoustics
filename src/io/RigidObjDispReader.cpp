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
    // first read displacement
    ReadDisplacement(fileDisplacement); 

    // read velocity 
    {
        std::cout << "Read velocity from rigidsim data : " << fileVelocity << "\n"; 
        std::ifstream inFile(fileVelocity, std::ios::in | std::ios::binary); 
        if (!inFile) throw std::runtime_error("**ERROR** Can't open file: " + fileVelocity); 

        REAL time; 
        if (!inFile.eof())
            inFile.read((char*)&time, sizeof(REAL));
        else 
            throw std::runtime_error("**WARNING** No time step found in the file"); 

        int v_totalObjects = 0; 
        int v_totalFrames = 0;
        while(true)
        {
            int objID = -1; 
            int totalObjects = 0; 
            VecVector3 frameVelocity; 
            VecVector3 frameAngularVelocity; 
            inFile.read((char*)&objID, sizeof(int)); 
            while(objID != -1) // end flag for the time step
            {
                Vector3<REAL> velocity; 
                Vector3<REAL> angularVelocity; 
                inFile.read((char*)&velocity, sizeof(Vector3<REAL>)); 
                inFile.read((char*)&angularVelocity, sizeof(Vector3<REAL>)); 
                frameVelocity.push_back(velocity); 
                frameAngularVelocity.push_back(angularVelocity); 
                totalObjects ++;
                inFile.read((char*)&objID, sizeof(int)); 
            }
            REAL timeCached = _timesteps.at(v_totalFrames); 
            if (fabs(timeCached - time) > 1E-10)
                throw std::runtime_error("**ERROR** velocity has different time frame than displacement."+std::to_string(timeCached)+"->"+std::to_string(time)); 

            _timesteps.push_back(time); 
            _velocity.push_back(frameVelocity); 
            _angularVelocity.push_back(frameAngularVelocity); 
            v_totalFrames ++; 
            if (totalObjects != 0) v_totalObjects = totalObjects; 

            if (!inFile.eof())
                inFile.read((char*)&time, sizeof(REAL));
            else 
                break; 
        }
        std::cout << " Frames read : " << v_totalFrames << std::endl;
        std::cout << " Objects read: " << v_totalObjects << std::endl;

        if (v_totalFrames != _totalFrames) 
            throw std::runtime_error("**ERROR** velocity has different number of frames than displacement"); 
        if (v_totalObjects != _totalObjects) 
            throw std::runtime_error("**ERROR** velocity has different number of frames than displacement"); 
    }
      
    // read acceleration
    {
        std::cout << "Read acceleration from rigidsim data : " << fileAcceleration << "\n"; 
        std::ifstream inFile(fileAcceleration, std::ios::in | std::ios::binary); 
        if (!inFile) throw std::runtime_error("**ERROR** Can't open file: " + fileAcceleration); 

        REAL time; 
        if (!inFile.eof())
            inFile.read((char*)&time, sizeof(REAL));
        else 
            throw std::runtime_error("**WARNING** No time step found in the file"); 

        int v_totalObjects = 0; 
        int v_totalFrames = 0;
        while(true)
        {
            int objID = -1; 
            int totalObjects = 0; 
            VecVector3 frameAcceleration; 
            VecVector3 frameAngularAcceleration; 
            inFile.read((char*)&objID, sizeof(int)); 
            while(objID != -1) // end flag for the time step
            {
                Vector3<REAL> acceleration; 
                Vector3<REAL> angularAcceleration; 
                inFile.read((char*)&acceleration, sizeof(Vector3<REAL>)); 
                inFile.read((char*)&angularAcceleration, sizeof(Vector3<REAL>)); 
                frameAcceleration.push_back(acceleration); 
                frameAngularAcceleration.push_back(angularAcceleration); 
                totalObjects ++;

                std::cout << objID << " : " << acceleration << ", " << angularAcceleration << std::endl;
                                                                    
                inFile.read((char*)&objID, sizeof(int)); 
            }
            REAL timeCached = _timesteps.at(v_totalFrames); 
            if (fabs(timeCached - time) > 1E-10)
                throw std::runtime_error("**ERROR** acceleration has different time frame than displacement."+std::to_string(timeCached)+"->"+std::to_string(time)); 

            _timesteps.push_back(time); 
            _acceleration.push_back(frameAcceleration); 
            _angularAcceleration.push_back(frameAngularAcceleration); 
            v_totalFrames ++; 
            if (totalObjects != 0) v_totalObjects = totalObjects; 

            if (!inFile.eof())
                inFile.read((char*)&time, sizeof(REAL));
            else 
                break; 
        }
        std::cout << " Frames read : " << v_totalFrames << std::endl;
        std::cout << " Objects read: " << v_totalObjects << std::endl;

        if (v_totalFrames != _totalFrames) 
            throw std::runtime_error("**ERROR** acceleration has different number of frames than displacement"); 
        if (v_totalObjects != _totalObjects) 
            throw std::runtime_error("**ERROR** acceleration has different number of frames than displacement"); 
    }

    std::cout << "Frames and objects check completed" << std::endl;
}
