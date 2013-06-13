#ifndef LINEARALGEBRA_PARDISO_SOLVER_HPP
#   define LINEARALGEBRA_PARDISO_SOLVER_HPP

#include <vector>
#include "Vector3.hpp"
#include "PardisoMatrix.hpp"

class PardisoSolver
{
    public:
        PardisoSolver():m_maxFact(1), m_num(1), 
                m_msgLvl(0), m_doSymbFact(true)
        {
            memset(m_parm, 0, sizeof(m_parm));
            memset(m_pt, 0, sizeof(m_pt));
        }

        /*! max number of factorization maintained by PARDISO */
        int& max_num_fact() { return m_maxFact; }

        int& message_level() { return m_msgLvl; }
        int* parameters() { return m_parm; }
        void reset_symbolic_analysis() { m_doSymbFact = true; }

        void solve(const PardisoMatrix<double>& A, 
                   const double *b, double *x);
        void release_memory();

    private:
        int     m_maxFact;
        int     m_num;
        int     m_type;
        int     m_msgLvl;       // message level
        int     m_parm[64];     // parameters to control the PARDISO solver
        int     m_size;         // size of the linear system
        void*   m_pt[64];       // memory pointer

        /*
         * Indicate if this solver is gonna do analysis and symbolic factorization
         * This a symbolic factorization is done, this variable becomes false. It 
         * never does symbolic factorization again until user explicitly set it (by
         * reset_symbolic_analysis)
         */
        bool    m_doSymbFact;   
};

#endif
