/******************************************************************************
 *  File: stvk.hpp
 *  Implement the St. Venant-Kirchhoff model
 *  Copyright (c) 2007 by Changxi Zheng
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#ifndef DEFORMABLE_STVK_H
#   define DEFORMABLE_STVK_H

#include "config.h"
#include "geometry/Tet.hpp"
#include "linearalgebra/Matrix3.hpp"
#include "linearalgebra/Matrix.hpp"

class StVKMaterial
{
    public:
        StVKMaterial() { }
        StVKMaterial(REAL lam, REAL mu):m_lambda(lam), m_mu(mu)
        { }
        StVKMaterial(const StVKMaterial& mat):
                m_lambda(mat.m_lambda), m_mu(mat.m_mu)
        { }
        StVKMaterial(const StVKMaterial* mat):
                m_lambda(mat->m_lambda), m_mu(mat->m_mu)
        { }

        void set_parameters(REAL l, REAL m)
        {
            m_lambda = l;
            m_mu = m;
        }

        /*! first Piola Kirchhoff stress tensor */
        void first_piola_kirchhoff(const Matrix3<REAL>& F, 
                Matrix3<REAL>& out) const;
        /*! second Piola Kirchhoff stress tensor */
        void second_piola_kirchhoff(const Matrix3<REAL>& F, 
                Matrix3<REAL>& out) const;      // $$TESTED
        void cauchy_stress_tensor(const Matrix3<REAL>& F,
                Matrix3<REAL>& out) const;
        /*! evaluate the stiffness matrix for a single tetrahedron */
        void stiffness_matrix(const Tet<REAL>& tet,
                Matrix<REAL>& out) const;       // $$TESTED
        /*! strain energy density in undeformed volumn */
        REAL strain_energy_density(const Matrix3<REAL>& F) const;

    protected:
        REAL    m_lambda;
        REAL    m_mu;
};

#endif
