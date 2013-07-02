#ifndef FRACTURE_SND_CONSTANTS_H
#   define FRACTURE_SND_CONSTANTS_H

#include "config.h"

const int COLLISION_DETECT_TOT_ITS  = 10;   // # of iterations for collision detection
const int COLLISION_DETECT_PAIR_ITS = 5;    // # of iterations for each pair of interpenetrated objects

const int CONTACT_DETECT_TOT_ITS    = 10;   // # of iterations for collision detection
const int CONTACT_DETECT_PAIR_ITS   = 20;   // # of iterations for each pair of interpenetrated objects
const int CONTACT_SHOCK_PROP_ITS    = 6;

const REAL CONT_REST_COEFF[] = { 
            -0.8, -0.4, -0.2, 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0., 0., 0., 0., 0., 0.,
              0.,   0.,   0., 0., 0.};

inline REAL contact_rest_coeff(int id)
{
    return id < 100 ? CONT_REST_COEFF[id] : 0.;
}

/*
 * The resolution of the level-set grids
 */
const REAL LEVELSET_RES = 80; //128 - 9;

const REAL GRAVITY      = 9.8;
#endif

