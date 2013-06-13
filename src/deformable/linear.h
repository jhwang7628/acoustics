#ifndef DEFORMABLE_LINEAR_H
#   define DEFORMABLE_LINEAR_H

#include "config.h"
#include "geometry/Tet.hpp"
#include "linearalgebra/Matrix.hpp"

class LinearMaterial
{
    public:
        LinearMaterial() { }
        LinearMaterial(REAL lam, REAL mu):m_lambda(lam), m_mu(mu)
        { }

        void set_parameters(REAL l, REAL m)
        {
            m_lambda = l;
            m_mu = m;
        }

        /*! 
         * evaluate the stiffness matrix for a single tetrahedron 
         * AT REST POSITION
         * Use the formula given in paper Brien et al. 2002 (Synthesizing Sounds from Rigid-Body Simulations)
         */
        void rest_stiffness_matrix(const Tet<REAL>& tet,
                Matrix<REAL>& out) const;
    private:
        REAL    m_lambda;
        REAL    m_mu;
};

#endif
