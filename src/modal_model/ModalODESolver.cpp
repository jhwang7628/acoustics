#include <iostream>
#include <cassert>
#include <modal_model/ModalODESolver.h> 

//##############################################################################
//##############################################################################
void ModalODESolver::
Initialize(ModalMaterialPtr &material, const REAL &omegaSquared, const REAL &timeStepSize)
{
    _material = material; 
    _omegaSquared = omegaSquared; 
    _omega = sqrt(omegaSquared);
    _timeStepSize = timeStepSize; 

    const REAL xi = material->xi(_omega); 
    if (xi > 1 || xi < 0) 
        throw std::runtime_error("**ERROR** xi is out of range [0,1]. check material parameters");

    const REAL omega_di = material->omega_di(_omega); 
    _epsilon = exp(-xi * _omega * timeStepSize); 
    _theta = omega_di * timeStepSize; 
    _gamma = asin(xi); 

    _2_epsilon_cosTheta = 2.0 * _epsilon * cos(_theta); 
    _epsilon_squared = pow(_epsilon, 2); 
    _coeff_Q_i = 2.0/(3.0*_omega*omega_di)
               * (_epsilon*cos(_theta+_gamma) - _epsilon_squared*cos(2.0*_theta+_gamma))
               * _material->inverseDensity; 

    _initialized = true; 
}

//##############################################################################
// Advance the state from (k-1) to (k). 
// Input:
//  qOld is q^(k-2) for this mode
//  qNew is q^(k-1) for this mode
// Output: 
//  qOld is q^(k-1) for this mode
//  qNew is q^(k)   for this mode
//##############################################################################
void ModalODESolver::
StepSystem(REAL &qOld, REAL &qNew, const REAL &Q)
{
    assert(_initialized); 
    const REAL q = _2_epsilon_cosTheta * qNew
                 - _epsilon_squared    * qOld
                 + _coeff_Q_i          * Q; 
    qOld = qNew; 
    qNew = q; 
    _time += _timeStepSize; 
}