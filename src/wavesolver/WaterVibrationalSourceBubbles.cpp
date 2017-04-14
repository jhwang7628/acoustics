#include <wavesolver/WaterVibrationalSourceBubbles.h>
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <sndgen/WavReader.hpp>
#include <utils/STL_Wrapper.h>
#include <fstream>

#include "bubbles/FileInput.hpp"

//##############################################################################
//##############################################################################
WaterVibrationalSourceBubbles::
WaterVibrationalSourceBubbles(RigidObjectPtr owner, const std::string &dataDir)
    : VibrationalSource(owner), _surfaceMesh(owner->GetMeshPtr())
{
    Initialize(dataDir);
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSourceBubbles::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle)
{
    if (normal.dotProduct(_wantedNormal)/normal.length() > _validAngleThreshold)
    {
        // TODO: return acceleration at this point
    }
    else
        return 0.0;
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSourceBubbles::
Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time)
{
    throw std::runtime_error("**ERROR** Cannot sample water vibrational source using vertexID");
    return 0.0;
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSourceBubbles::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented");
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSourceBubbles::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented");
}

//##############################################################################
//##############################################################################
void WaterVibrationalSourceBubbles::
Initialize(const std::string &dataDir)
{
    std::cout << "Initialize WaterVibrationalSourceBubbles from directory: " << dataDir << std::endl;
    parseFileNames(dataDir);
}

