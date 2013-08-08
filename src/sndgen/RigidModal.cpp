#include "RigidModal.h"
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <string>
#include <fstream>
#include "SndGenConfig.h"
#include "utils/print_msg.h"

using namespace std;

RigidModal::RigidModal(libconfig::Setting& s)
{
    using namespace libconfig;

    string txt;
    if ( !s.lookupValue("eigenmodes", txt) )
    {
        PRINT_ERROR("Cannot find the configuration 'eigenmodes'\n");
        exit(3);
    }

    //// load modal paramters
    if ( !s.lookupValue("density", m_density) || 
         !s.lookupValue("alpha", m_alpha)     ||
         !s.lookupValue("beta",  m_beta) )
    {
        PRINT_ERROR("Cannot find 'density', 'alpha' or 'beta'\n");
        exit(3);
    }
    if ( !s.lookupValue("dampingGamma", m_gamma) )
    {
      PRINT_ERROR( "No inverse frequency damping 'dampingGamma' specified\n" );
      m_gamma = 0.0;
    }
    if ( !s.lookupValue("id", m_id) )
    {
      PRINT_ERROR( "Cannot find 'id'\n" );
      exit(3);
    }
    printf("Material density=%lf, alpha=%lf, beta=%lf, gamma = %lf\n",
            m_density, m_alpha, m_beta, m_gamma);
    m_invDensity = 1. / m_density;

    printf("Load eigenmodes [%s] ...\n", txt.c_str());
    load_eigenmodes(txt.c_str());
}

/*!
 * Load eigen modes of rigid object from file 
 */
void RigidModal::load_eigenmodes(const char* file)
{
    ifstream fin(file, ios::binary);

    fin.read((char *)&m_n3, sizeof(int));       // size of eigen problem
    fin.read((char *)&m_numModes, sizeof(int)); // number of modes
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot read file: %s\n", file);
        exit(4);
    }
    printf("Load eigen-modes %d, %d\n", m_numModes, m_n3);

    m_eigenmodes.resize(m_numModes);
    fin.read((char*)&m_eigenmodes[0], sizeof(double)*m_numModes);
    printf("Compute eigen-modes' frequencies ...\n");
    int nmds = 0;
    for(nmds = 0;nmds < m_numModes;++ nmds)
    {
        m_eigenmodes[nmds] *= m_invDensity;
        if ( m_eigenmodes[nmds] > CUTTING_OMEGA*CUTTING_OMEGA ) break;
    }
    printf("%d modes in audible range\n", nmds);
    m_numModes = nmds;

    //// all the eigen vectors are stored in a n3 x nModes matrix
    //// it is stored as v1 v2 ...
    m_eigenvec.resize(m_n3*m_numModes);
    fin.read((char *)&m_eigenvec[0], sizeof(double)*m_eigenvec.size());
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot read file: %s\n", file);
        exit(4);
    }
    fin.close();

    m_omega.resize(m_numModes);
    m_omegad.resize(m_numModes);
    m_xi.resize(m_numModes);
    m_wides.resize(m_numModes);
    for(int i = 0;i < m_numModes;++ i)
    {
        m_omega[i] = sqrt(m_eigenmodes[i]);
        if ( m_omega[i] < 125.6637 ) continue;

        // FIXME: This doesn't match equatino (10) in the IIR paper!
        double c = m_alpha * m_eigenmodes[i] + m_beta;
        m_xi[i] = c / (2. * m_omega[i]);
#if 0
        double c = m_alpha + m_beta * m_eigenmodes[i];
        m_xi[i] = c / ( 2.0 * m_omega[i] );
#endif
        if ( m_xi[i] >= 1. || m_xi[i] < 0. )
        {
            PRINT_ERROR("xi is supposed to be in the range [0, 1]: %lf\n", m_xi[i]);
            //exit(5);
        }
        m_omegad[i] = m_omega[i] * sqrt(1 - m_xi[i]*m_xi[i]);
        m_wides[i] = 2. * M_PI / m_omegad[i] * 0.5;
#if 0
        printf( "m_wides[ %d ] = %f\n", i, m_wides[ i ] );
#endif
    }
}

void RigidModal::modal_impulses(const int* impvtx, const double* imp, 
                                int impLen, double* modalImp) const
{   
#ifdef USE_OPENMP
    #pragma omp parallel for default(none) shared(modalImp, impLen, imp, impvtx)
#endif
    for(int im = 0;im < m_numModes;++ im)
    {
        const double* peigvec = &m_eigenvec[im * m_n3];
        modalImp[im] = 0;
        for(int i = 0, k = 0;i < impLen;++ i, k += 3)
        {
            const int iidd = impvtx[i]*3;
            modalImp[im] += imp[k]   * peigvec[iidd]  ;
            modalImp[im] += imp[k+1] * peigvec[iidd+1];
            modalImp[im] += imp[k+2] * peigvec[iidd+2];
        }
        modalImp[im] *= m_invDensity;
    }
}


