#include <wavesolver/ChirpVibrationalSource.h>

//##############################################################################
//##############################################################################
REAL ChirpVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle)
{
    if (time<_time[0])
        return 0.0;
    const REAL alpha = (time<_time[1] ?
            0.5 * (time - _time[0])/(_time[1] - _time[0]) :
            1.0);
    const REAL phase = (alpha < 0.52 ? 0 : (_time[1] ) * 0.5 * (_omega[1] - _omega[0]));
    const REAL omega =  (1.0 - alpha)*_omega[0] + alpha*_omega[1];
    //const REAL amp   = ((1.0 - alpha)*_amp[0]   + alpha*_amp[1])*pow(omega,2);
    const REAL amp = 1.0;
    //COUT_SDUMP(omega/2.0/M_PI);
    return sin(omega*time + phase)*amp;
}

//##############################################################################
//##############################################################################
REAL ChirpVibrationalSource::
Evaluate(const int &vertexID, const Vector3d &normal, const REAL &time)
{
    const REAL an = Evaluate(Vector3d(), normal, time);
    return an;
}

//##############################################################################
//##############################################################################
Vector3d ChirpVibrationalSource::
Evaluate(const int &vertexID, const REAL &time)
{
    Vector3d n = _owner->GetMeshPtr()->normal(vertexID);
    n = _owner->ObjectToWorldVector(n).normalized();
    const REAL an = Evaluate(vertexID, n, time);
    return n*an;
}
