#ifndef FAST_BEM_TRANSFER_INTERP_H
#   define FAST_BEM_TRANSFER_INTERP_H

#include <libconfig.h++>
#include <vector>
#include "interp/CSpline.hpp"
#include "utils/print_msg.h"

#include "linearalgebra/Vector3.hpp"

class FBemTransferInterp
{
    public:
        typedef CSpline<double, false>  TInterp;

        FBemTransferInterp(libconfig::Setting& s, int nmodes, 
                           const std::vector<double>& omega,
                           int objectID,
                           const Vector3d *centerOfMass = NULL);

        void init_interp(int mid);
        double transfer_norm(double ts);
        double transfer_distance(double ts);
        void stereo_transfer_norm(double, double&, double&)
        {
            PRINT_ERROR("FBemTransferInterp::stereo_transfer_norm is not supported!\n");
            exit(1);
        }

    private:
        void load_transfer(const char* fileptn, const std::vector<double>& omega,
                           int objectID);
        void load_time_ticks(const char* file, const Vector3d *centerOfMass);

        void save_transfer(const char* fileptn, int objectID);
        bool load_transfer_binary(const char* fileptn, int objectID);

    private:
        int     m_nModes;
        int     m_curMId;
        std::vector< double >               m_transferTs;
        std::vector< std::vector<double> >  m_transferNrm;  // norm of the transfer function value

        // Distances from the center of mass at each tick
        std::vector< double >               m_distances;
        TInterp                             m_distanceInterp;

        TInterp                             m_transferInterp;
};

#endif
