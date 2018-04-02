#ifndef CHIRP_VIBRATIONAL_SOURCE_H
#define CHIRP_VIBRATIONAL_SOURCE_H

#include <TYPES.h>
#include <wavesolver/VibrationalSource.h>
#include <wavesolver/FDTD_RigidObject.h>

//##############################################################################
// Chirp vibrations.
//##############################################################################
class ChirpVibrationalSource : public VibrationalSource
{
    private:
        REAL _omega[2];
        REAL _time[2];
        REAL _amp[2];

    public:
        ChirpVibrationalSource(RigidObjectPtr owner,
                                  const REAL &omega_0, const REAL &omega_1,
                                  const REAL &t_0    , const REAL &t_1,
                                  const REAL &a_0    , const REAL &a_1)
            : VibrationalSource(owner),
              _omega{omega_0, omega_1},
              _time{t_0, t_1},
              _amp{a_0, a_1}
        {}

        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle=-1);
        virtual REAL Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time);
        virtual Vector3d Evaluate(const int &vertexID, const REAL &time);
};

#endif
