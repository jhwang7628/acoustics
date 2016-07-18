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
//##############################################################################
void ModalODESolver::
StepSystem(const REAL &Q)
{
    assert(_initialized); 
    const REAL q = _2_epsilon_cosTheta * _qNew
                 - _epsilon_squared    * _qOld
                 + _coeff_Q_i          * Q; 
    _qOld = _qNew; 
    _qNew = q; 
}
