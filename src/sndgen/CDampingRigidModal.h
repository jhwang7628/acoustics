#ifndef CDAMPING_RIGID_MODEL_H
#   define CDAMPING_RIGID_MODEL_H

#include "RigidModal.h"

class CDampingRigidModal : public RigidModal
{
    public:
        CDampingRigidModal(libconfig::Setting& s) : RigidModal(s)
        {
            using namespace libconfig;
            
            if ( !s.lookupValue("gamma", gamma_) )
            {
                PRINT_ERROR("Cannot find the configuration 'gamma'\n");
                exit(3);
            }

            update_damping_coeff();
        }

        const std::vector<double>& cd_omegad() const
        {   return cdOmegad_; }
        const std::vector<double>& cd_xi() const
        {   return cdXi_; }

    private:
        // compute the damping coefficients with contact damping 
        void update_damping_coeff()
        {
            cdOmegad_.resize(m_numModes);
            cdXi_.resize(m_numModes);
            for(int i = 0;i < m_numModes;++ i)
            {
                if ( m_omega[i] < 125.6637 ) continue;

                double c = m_alpha * m_eigenmodes[i] + m_beta + gamma_;
                cdXi_[i] = c / (2. * m_omega[i]);
                if ( cdXi_[i] >= 1. || cdXi_[i] < 0. )
                {
                    PRINT_ERROR("xi_cd is supposed to be in the range [0, 1]: %lf\n", cdXi_[i]);
                    exit(5);
                }
                cdOmegad_[i] = m_omega[i] * sqrt(1. - cdXi_[i]*cdXi_[i]);
            }
        }

    private:
        double              gamma_;
        std::vector<double> cdOmegad_;
        std::vector<double> cdXi_;
};

#endif
