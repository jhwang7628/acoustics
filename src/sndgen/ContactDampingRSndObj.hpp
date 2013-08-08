#ifndef CONTACT_DAMPING_RIGID_SOUND_OBJ_HPP
#   define CONTACT_DAMPING_RIGID_SOUND_OBJ_HPP

#include "config.h"
#include <vector>
#include "RigidSoundObj.hpp"
#include "filter/LinearFilter.hpp"

/*
 * This class is for generating modal sound from rigid body simulator.
 * The contact damping is considered when synthesizing the sound.
 *
 * TModalAnalysis: modal analysis
 *              int num_modes()
 *              const vector<double>& omega()
 *              const vector<double>& omegad()
 *              const vector<double>& xi()
 * TImpulseEval:   evaluate modal impulses at given time
 *              has_impulses()
 * TTransferEval:  evaluate transfer function values for given frequency at given time
 */
template<class TModalAnalysis, 
         class TImpulseEval, 
         class TTransferEval>
class ContactDampingRSndObj : public RigidSoundObj<TModalAnalysis, TImpulseEval, TTransferEval>
{
    public:
        typedef SoundBuffer             TSndBuf;
        typedef StereoSndBuffer         TStereoBuf;
        typedef RigidSoundObj<TModalAnalysis,TImpulseEval,TTransferEval> TSndObj;

        ContactDampingRSndObj(libconfig::Setting& s, TModalAnalysis* pmodel,
                              TImpulseEval* pimp, TTransferEval* ptrans,
                              double d = 0);
        /*
         * generate audio signals for this object
         */
        int generate_sound(TSndBuf& buf); 

    private:
        typedef QuarticSplineKernel<double>     TCDFilter;

        int         cdThresh_;
        double      cdWinTimeLen_;
        TCDFilter   cdFilter_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////

template<class TModalAnalysis, 
         class TImpulseEval, 
         class TTransferEval>
ContactDampingRSndObj<TModalAnalysis, TImpulseEval, TTransferEval>::ContactDampingRSndObj(
        libconfig::Setting& s,
        TModalAnalysis* pmodel, TImpulseEval* pimp, TTransferEval* ptrans,
        double d): 
        RigidSoundObj<TModalAnalysis, TImpulseEval, TTransferEval>(s, pmodel, pimp, ptrans, d)
{
    using namespace libconfig;

    int winlen;
    if ( !s.lookupValue("contact_damping_winlen", winlen) )
    {
        PRINT_ERROR("Cannot find the variable path for 'contact_damping_winlen'\n");
        exit(3);
    }

    if ( !s.lookupValue("contact_damping_threshold", cdThresh_) )
    {
        PRINT_ERROR("Cannot find the variable path for 'contact_damping_threshold'\n");
        exit(3);
    }

    cdWinTimeLen_ = (double)winlen * pimp->sim_step_size();
    cdFilter_.set_win_len(cdWinTimeLen_*2.);
}

template<class TModalAnalysis, 
         class TImpulseEval, 
         class TTransferEval>
int ContactDampingRSndObj<TModalAnalysis, TImpulseEval, TTransferEval>::generate_sound(TSndBuf& buf)
{
    using namespace std;

    if ( !this->mp_imp->has_impulses() ) return SUCC_RETURN;

    const double H = 1. / (double)buf.sample_rate();
    const int NUM_MODES = this->mp_modal->num_modes();
    const vector<double>& omega  = this->mp_modal->omega();
    const vector<double>& omegad = this->mp_modal->omegad();
    const vector<double>& xi     = this->mp_modal->xi();
    const vector<double>& gwide  = this->mp_modal->gaussian_wide();
    const vector<double>& cd_omegad = this->mp_modal->cd_omegad();
    const vector<double>& cd_xi  = this->mp_modal->cd_xi();

    for(int mid = 0;mid < NUM_MODES;++ mid)
    {
        if ( omega[mid] < 125.6637 || omega[mid] > CUTTING_OMEGA ) continue;
        printf("CDAMPING::freq[%d] = %lf\n", mid, omega[mid]*0.5*M_1_PI);

        // initialize transfer for this mode
        this->mp_transfer->init_interp(mid);

        const double WIDE    = fmin(gwide[mid], 25.*H);
        const double TSBEGIN = fmax(0., this->mp_imp->begin_time() - 4*WIDE);
        const double TSEND   = this->mp_imp->end_time() - log(1E-5)/(xi[mid] * omega[mid]);
        this->mp_imp->init_interp(mid, WIDE);

        int LL = 0;
        double lastm1 = 0, lastm2 = 0, cd_lastm1 = 0, cd_lastm2 = 0;
        //// IIR filter
        double theta = omegad[mid] * H;
        double cd_theta = cd_omegad[mid] * H;
        double eps    = exp(-xi[mid] * omega[mid] * H);
        double cd_eps = exp(-cd_xi[mid] * omega[mid] * H);
        double gamma = asin(xi[mid]);
        double cd_gamma = asin(cd_xi[mid]);
        double c1 = 2 * eps * cos(theta),
               c2 = - M_SQR(eps),
               c3 = 2*(eps*cos(theta+gamma) + c2*cos(2*theta+gamma)) / (3*omega[mid]*omegad[mid]);
        double cd_c1 = 2 * cd_eps * cos(cd_theta),
               cd_c2 = - M_SQR(cd_eps),
               cd_c3 = 2*(cd_eps*cos(cd_theta+cd_gamma) + cd_c2*cos(2*cd_theta+cd_gamma)) / (3*omega[mid]*cd_omegad[mid]);

        buf.begin_sound_gen(TSBEGIN);
        double ticktime;
        int vlow, vhigh;
        for(int tick = 0;(ticktime = TSBEGIN + tick*H) < TSEND;++ tick)
        {
            double I = this->mp_imp->get_impulse(ticktime);
            double u    = c1*lastm1 + c2*lastm2 + c3*I;
            double cd_u = cd_c1*cd_lastm1 + cd_c2*cd_lastm2 + cd_c3*I;

            this->mp_imp->current_interval_info(vlow, vhigh);
            LL = std::max(LL, vhigh-vlow);
            bool bCD = (vhigh - vlow) > cdThresh_ && ticktime > this->mp_imp->impulse_ts(vlow);
            if ( !bCD )
            {
                if ( !buf.add_sound_sample(this->mp_transfer->transfer_norm(ticktime)*u) ) break;
            }
            else
            {
                double timed = ticktime - this->mp_imp->impulse_ts(vlow);
                if ( timed >= cdWinTimeLen_ )
                {
                    if ( !buf.add_sound_sample(this->mp_transfer->transfer_norm(ticktime)*cd_u) ) break;
                }
                else
                {
                    double w = cdFilter_(cdWinTimeLen_ - timed);
                    if ( !buf.add_sound_sample(
                                this->mp_transfer->transfer_norm(ticktime) *
                                (cd_u*w + (1.-w)*u)) )
                        break;
                }
            }
            lastm2    = lastm1;     lastm1    = u;
            cd_lastm2 = cd_lastm1;  cd_lastm1 = cd_u;
        }
        printf("LL = %d\n", LL);
    }
    return SUCC_RETURN;
}

#endif
