#ifndef WATER_VIBRATIONAL_SOURCE_BUBBLES_H
#define WATER_VIBRATIONAL_SOURCE_BUBBLES_H

#include <Eigen/Dense>
#include <TYPES.h>
#include <wavesolver/VibrationalSource.h>
#include <wavesolver/FDTD_RigidObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>

#include "bubbles/Oscillator.hpp"
#include "bubbles/Mesh.hpp"

//##############################################################################
// This class handles source evaluation for water surface objects using data from the bubbles project.
//##############################################################################
class WaterVibrationalSource : public VibrationalSource
{
    public:
        typedef std::shared_ptr<TriangleMesh<REAL> > TriangleMeshPtr;

    private:
        Vector3d            _wantedNormal = Vector3d(0, 1, 0);
        REAL                _validAngleThreshold = 0.5; // use to determine if the vertex has source. See Evaluate() for usage. default: 0.5
        TriangleMeshPtr     _surfaceMesh;
        REAL                _sampleRate;
        REAL                _startTime = 0.0;

        REAL                _dt;

    public:
        WaterVibrationalSourceBubbles(RigidObjectPtr owner, const std::string &dataDir);

        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle);
        virtual REAL Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time);
        virtual REAL EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time);
        virtual REAL EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time);

        void Initialize(const std::string &dataDir);
};

#endif
