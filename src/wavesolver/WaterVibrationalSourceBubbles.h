#ifndef WATER_VIBRATIONAL_SOURCE_BUBBLES_H
#define WATER_VIBRATIONAL_SOURCE_BUBBLES_H

#include <Eigen/Dense>
#include <TYPES.h>
#include <wavesolver/VibrationalSource.h>
#include <wavesolver/FDTD_RigidObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <geometry/KDTree.hpp>

#include "bubbles/Bubble.hpp"
#include "bubbles/Oscillator.hpp"
#include "bubbles/Mesh.hpp"
#include "bubbles/FileInput.hpp"

//##############################################################################
// This class handles source evaluation for water surface objects using data from the bubbles project.
//##############################################################################
class WaterVibrationalSourceBubbles : public VibrationalSource
{
    public:
        typedef std::shared_ptr<TriangleMesh<REAL> > TriangleMeshPtr;

        // indexed by bubble id, then vector of triangle velocities
        typedef MLSModeInterpolator<double, 3, 1> MLSInterp; // TODO: should this be 1d or 3d? (interpolate normal velocities or full velocity vectors?)

        struct DistSq
        {
            double operator() (const Eigen::Vector3d &a, const Eigen::Vector3d &b) {return (a-b).squaredNorm();}
        };

        typedef KDTree<3, Eigen::Vector3d, DistSq> PointKDTree;

        enum FreqType
        {
            CAPACITANCE,
            MINNAERT
        };

    private:
        Vector3d            _wantedNormal = Vector3d(0, 1, 0);
        REAL                _validAngleThreshold = 0.5; // use to determine if the vertex has source. See Evaluate() for usage. default: 0.5
        TriangleMeshPtr     _surfaceMesh;

        REAL                _dt;

        double _curTime; // current time step
        double _t1, _t2; // surrounding times for surface data
        Mesh _m1, _m2;
        Mesh *_m; // mesh that is currently being used
        SurfaceVelocityData _v1, _v2;
        MLSInterp _mls;
        std::shared_ptr<PointKDTree> _kd1, _kd2; // TODO: add copy/move semantics to the kd tree class so shared pointers aren't necessary
        std::vector<BubbleInputInfo> _b1, _b2;
        Eigen::VectorXd _velT1, _velT2; // the total velocities before and after current time (not at _t1 and _t2)

        Eigen::VectorXd _accel;
        Eigen::VectorXd _projectedAccel;

        std::map<double, FileNames> _fileInfo; // indexed by time
        std::map<int, Bubble> _bubbles;
        std::vector<Oscillator> _oscillators;

        void step(REAL time);
        void parseConfigFile(const std::string &, FreqType fType);
        std::pair<int, Bubble> parseBubbleInfo (std::ifstream &in, FreqType fType);

        // Attach the bubbles together to get the full oscillator lifetime
        void makeOscillators(const std::map<int, Bubble> &singleBubbles);

        void updateOscillators(REAL time);

        void computeVelocities(REAL time);

        void projectToSurface();

    public:
        WaterVibrationalSourceBubbles(RigidObjectPtr owner, const std::string &dataDir);

        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle);
        virtual REAL Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time);
        virtual REAL EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time);
        virtual REAL EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time);

        void Initialize(const std::string &dataDir);
};

#endif
