/*
 * =====================================================================================
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
 * 
 * -------------------------------------------------------------------------------------
 *
 *       Filename:  JacobiEigenSolve3x3.hpp
 *
 *    Description:  Eigen solve of 3x3 matrix using Jacobi iteration
 *
 *        Version:  1.0
 *        Created:  09/02/11 00:44:29
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#ifndef JACOBI_EIGEN_SOLVE_3X3_INC
#   define JACOBI_EIGEN_SOLVE_3X3_INC

template <typename T>
class JacobiEigenSolve3x3
{
    public:
        JacobiEigenSolve3x3();

        /*! Eigen solve */
        void solve(const T A[3][3], T V[3][3], T d[3]);

    private:
};

#endif
