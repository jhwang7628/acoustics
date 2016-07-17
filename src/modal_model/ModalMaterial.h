#ifndef MODAL_MATERIAL_H 
#define MODAL_MATERIAL_H 
#include <config.h>

//##############################################################################
// This struct stores and computes all the material paramters needed for modal
// analysis. The naming convention follows DyRT paper [James 2002]. In
// particular, IIR filter equation is useful to time-step the modal equations.
// (eq. 11, 19).
//
// Can be built by parser/ImpulseResponseParser. 
//##############################################################################
struct ModalMaterial
{
    int id; 
    REAL alpha; 
    REAL beta; 
    REAL density; // assuming we always have constant, uniform density
    REAL inverseDensity; // cached
    inline REAL xi(const REAL &omega_i){return 0.5*(alpha/omega_i + beta*omega_i);} // eq.10, xi = [0,1]
    inline REAL omega_di(const REAL &omega_i){return omega_i*sqrt(1.0 - pow(xi(omega_i),2));} // eq.12.
};

#endif
