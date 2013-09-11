#ifndef SOUND_RIGID_SND_OBJ_H
#   define SOUND_RIGID_SND_OBJ_H

#include "SndGenConfig.h"
#include <vector>
#ifdef USE_MKL
#   include <mkl.h>
#else
#   if defined(__APPLE__) || defined(MACOSX)
#       include <vecLib/cblas.h>
#   endif
#endif
#include <libconfig.h++>

#include <linearalgebra/Vector3.hpp>

#include "SoundBuffer.hpp"

#include <utils/math.hpp>
#include <utils/STLUtil.h>

/*
 * This class is for generating modal sound from rigid body simulator.
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
    class RigidSoundObj
{
    public:
        typedef SoundBuffer             TSndBuf;
        typedef StereoSndBuffer         TStereoBuf;

        /*!
         * \param d     sound delay
         */
        RigidSoundObj(libconfig::Setting& s, TModalAnalysis* pmodel, 
                TImpulseEval* pimp, TTransferEval* ptrans,
                double d = 0);

        // If we want to use a hertz impulse series
        RigidSoundObj(libconfig::Setting& s, TModalAnalysis* pmodel, 
                TImpulseEval* pimp,
                HertzImpulseSeries<TModalAnalysis> *hertzImpulseSeries,
                TTransferEval* ptrans,
                double d = 0);

        //virtual ~RigidSoundObj();

        /*
         * generate audio signals for this object
         */
        virtual int generate_sound(TSndBuf& buf); 
        virtual int generate_hertz_sound(TSndBuf& buf);
        virtual int timestep_hertz_sound(TSndBuf& buf, bool useDelay = false);
        virtual int generate_stereo_snd(TStereoBuf&);

        int obj_id() const
        {   return m_objId; }

        /* load impulses */
        void begin_read_impulse()
        {   mp_imp->begin_read_impulse(mp_modal); }
        void add_impulse(double ts, int vtxid, const Vector3d& imp, char onsurf)
        {   mp_imp->add_impulse(ts, vtxid, imp, onsurf); }
        void add_impulse(double ts, int vtxid, const Vector3d &imp,
                double impulseLength, double impulseScale)
        {
            mp_hertzImpulseSeries->add_impulse( vtxid, ts,
                    impulseLength, impulseScale,
                    imp );
        }
        void end_read_impulse()
        {   mp_imp->end_read_impulse(); }

        bool load_binary_impulses( const char *filename )
        {
            return mp_hertzImpulseSeries->load( filename );
        }

        TModalAnalysis* modal_analysis() 
        {   return mp_modal; }
        TImpulseEval*   impulse_interp() 
        {   return mp_imp; }
        TTransferEval*  transfer_interp() 
        {   return mp_transfer; }
        HertzImpulseSeries<TModalAnalysis> *hertz_series()
        {   return mp_hertzImpulseSeries; }

    protected:
        int         m_objId;
        double      m_delay;

        TModalAnalysis* mp_modal;
        TImpulseEval*   mp_imp; 
        TTransferEval*  mp_transfer;

        HertzImpulseSeries<TModalAnalysis> *mp_hertzImpulseSeries;
};

/////////////////////////////////////////////////////////////////////
template<class TModalAnalysis, 
    class TImpulseEval, 
    class TTransferEval>
RigidSoundObj<TModalAnalysis, TImpulseEval, TTransferEval>::RigidSoundObj(
        libconfig::Setting& s, 
        TModalAnalysis* pmodel, TImpulseEval* pimp, TTransferEval* ptrans,
        double d)
: m_delay(d),
    mp_modal(pmodel),
    mp_imp(pimp),
    mp_transfer(ptrans),
    mp_hertzImpulseSeries( NULL )
{
    using namespace libconfig;
    if ( !s.lookupValue("id", m_objId) )
    {
        PRINT_ERROR("Cannot find property 'id'\n");
        exit(3);
    }
}

template<class TModalAnalysis, 
    class TImpulseEval, 
    class TTransferEval>
RigidSoundObj<TModalAnalysis, TImpulseEval, TTransferEval>::RigidSoundObj(
        libconfig::Setting& s, 
        TModalAnalysis* pmodel,
        TImpulseEval* pimp,
        HertzImpulseSeries<TModalAnalysis> *hertzImpulseSeries,
        TTransferEval* ptrans,
        double d)
: m_delay(d),
    mp_modal(pmodel),
    mp_imp(pimp),
    mp_transfer(ptrans),
    mp_hertzImpulseSeries( hertzImpulseSeries )
{
    using namespace libconfig;
    if ( !s.lookupValue("id", m_objId) )
    {
        PRINT_ERROR("Cannot find property 'id'\n");
        exit(3);
    }
}

template<class TModalAnalysis, 
    class TImpulseEval, 
    class TTransferEval>
    int RigidSoundObj<TModalAnalysis, TImpulseEval,
    TTransferEval>::generate_sound(TSndBuf& buf)
{
    using namespace std;

    if ( !mp_imp->has_impulses() ) return SUCC_RETURN;

    const double H = 1. / (double)buf.sample_rate();
    const int NUM_MODES = mp_modal->num_modes();
    const vector<double>& omega  = mp_modal->omega();
    const vector<double>& omegad = mp_modal->omegad();
    const vector<double>& xi     = mp_modal->xi();
    const vector<double>& gwide  = mp_modal->gaussian_wide();

    for(int mid = 0;mid < NUM_MODES;++ mid)
    {
        if ( omega[mid] < 125.6637 || omega[mid] > CUTTING_OMEGA ) continue;
        printf("freq[%d] = %lf\n", mid, omega[mid]*0.5*M_1_PI);

        // initialize transfer for this mode
        mp_transfer->init_interp(mid);

        const double WIDE    = fmin(gwide[mid], 25.*H);
        const double TSBEGIN = fmax(0., mp_imp->begin_time() - 4*WIDE);
        const double TSEND   = mp_imp->end_time()
            - log(1E-5)/(xi[mid] * omega[mid]);
        mp_imp->init_interp(mid, WIDE);

        double lastm1 = 0, lastm2 = 0;
        //// IIR filter
        double theta = omegad[mid] * H;
        double eps   = exp(-xi[mid] * omega[mid] * H);
        double gamma = asin(xi[mid]);
        double c1 = 2 * eps * cos(theta),
               c2 = - M_SQR(eps),
               c3 = 2*(eps*cos(theta+gamma) + c2*cos(2*theta+gamma))
                   / (3*omega[mid]*omegad[mid]);

        buf.begin_sound_gen(TSBEGIN);
        double ticktime;
        for(int tick = 0;(ticktime = TSBEGIN + tick*H) < TSEND;++ tick)
        {
            double u = c1*lastm1 + c2*lastm2 + c3*mp_imp->get_impulse(ticktime);
            if ( !buf.add_sound_sample(mp_transfer->transfer_norm(ticktime)*u) )
                break;
            lastm2 = lastm1;
            lastm1 = u;
        }
    }
    return SUCC_RETURN;
}

//////////////////////////////////////////////////////////////////////
// Generates sound using a hertz contact model
//////////////////////////////////////////////////////////////////////
template<class TModalAnalysis, 
    class TImpulseEval, 
    class TTransferEval>
    int RigidSoundObj<TModalAnalysis, TImpulseEval,
    TTransferEval>::generate_hertz_sound(TSndBuf& buf)
{
    using namespace std;

#if 0
    if ( !mp_imp->has_impulses() ) return SUCC_RETURN;
#endif

    const double H = 1. / (double)buf.sample_rate();
    const int NUM_MODES = mp_modal->num_modes();
    const vector<double>& omega  = mp_modal->omega();
    const vector<double>& omegad = mp_modal->omegad();
    const vector<double>& xi     = mp_modal->xi();

    for(int mid = 0;mid < NUM_MODES;++ mid)
    {
        if ( omega[mid] < 125.6637 || omega[mid] > CUTTING_OMEGA ) continue;
        printf("freq[%d] = %lf\n", mid, omega[mid]*0.5*M_1_PI);

        // initialize transfer for this mode
        mp_transfer->init_interp(mid);

        const double TSBEGIN = fmax( 0., mp_hertzImpulseSeries->startTime() );
        const double TSEND   = mp_hertzImpulseSeries->endTime();

        mp_hertzImpulseSeries->initInterp( mid );

        double lastm1 = 0, lastm2 = 0;
        //// IIR filter
        double theta = omegad[mid] * H;
        double eps   = exp(-xi[mid] * omega[mid] * H);
        double gamma = asin(xi[mid]);
        double c1 = 2 * eps * cos(theta),
               c2 = - M_SQR(eps),
               c3 = 2*(eps*cos(theta+gamma) + c2*cos(2*theta+gamma))
                   / (3*omega[mid]*omegad[mid]);

        buf.begin_sound_gen(TSBEGIN);
        double ticktime;
        for(int tick = 0;(ticktime = TSBEGIN + tick*H) < TSEND;++ tick)
        {
            //double u = c1*lastm1 + c2*lastm2 + c3*mp_imp->get_impulse(ticktime);
            double u = c1*lastm1 + c2*lastm2
                + c3*mp_hertzImpulseSeries->getForce(ticktime);
            if ( !buf.add_sound_sample(mp_transfer->transfer_norm(ticktime)*u) )
                break;
            lastm2 = lastm1;
            lastm1 = u;
        }
    }
    return SUCC_RETURN;
}

//////////////////////////////////////////////////////////////////////
// Generates sound using a hertz contact model and explicit Newmark
// time stepping
//////////////////////////////////////////////////////////////////////
template<class TModalAnalysis, 
    class TImpulseEval, 
    class TTransferEval>
    int RigidSoundObj<TModalAnalysis, TImpulseEval,
    TTransferEval>::timestep_hertz_sound(TSndBuf& buf,
            bool useDelay)
{
    using namespace std;

#if 0
    if ( !mp_imp->has_impulses() ) return SUCC_RETURN;
#endif

    const double H = 1. / (double)buf.sample_rate();
    const double H2 = H * H;
    const int NUM_MODES = mp_modal->num_modes();
    const vector<double>& omega  = mp_modal->omega();
    const vector<double>& omegad = mp_modal->omegad();
    const vector<double>& xi     = mp_modal->xi();

    double alpha = mp_modal->alpha();
    double beta = mp_modal->beta();
    double gamma = mp_modal->gamma();

    double delay;

#ifdef DEBUG
    FloatArray forceProfile;
#endif

    printf( "Generating sound from %d modes\n", NUM_MODES );

    for(int mid = 0;mid < NUM_MODES;++ mid)
    {
        printf("Time step freq[%d] = %lf\n", mid, omega[mid]*0.5*M_1_PI);

#ifdef DEBUG
        forceProfile.clear();
#endif

        // initialize transfer for this mode
        mp_transfer->init_interp( mid );
        mp_hertzImpulseSeries->initInterp( mid );

        double                 accelOld = 0.0;
        double                 accel = 0.0;
        double                 displacement = 0.0;
        double                 velocity = 0.0;
        double                 force;
        double                 vCoefficient;
        double                 newmarkVCoefficient;

        double                 omega2 = omega[ mid ] * omega[ mid ];
        //double                 damping = alpha + beta * omega2;
        // FIXME: Reversing the meaning of alpha and beta here, since
        // we have done this in all of our configuration file
        double                 damping = alpha * omega2 + beta;
        double                 f = omega[ mid ] / ( 2 * M_PI );

        //cout << "damping = " << damping << endl;
        //cout << "gamma = " << gamma << endl;

        // FIXME: Trying out some alternative damping here
        damping += gamma / omega2;
        //cout << "damping = " << damping << endl;

        vCoefficient = ( 1.0 + damping * H + omega2 * H * H );
        newmarkVCoefficient = ( 1.0 + damping * H / 2.0 );

        //cout << "vCoefficient = " << vCoefficient << endl;

        const double TSBEGIN = fmax( 0., mp_hertzImpulseSeries->startTime() );
        const double TSEND   = mp_hertzImpulseSeries->endTime();

        double ticktime;
        buf.begin_sound_gen(TSBEGIN);
        for(int tick = 0;(ticktime = TSBEGIN + tick*H) < TSEND;++ tick)
        {
            // Time delay
            delay = mp_transfer->transfer_distance( ticktime ) / 343.0;

            force = mp_hertzImpulseSeries->getForce( ticktime );

#ifdef DEBUG
            forceProfile.push_back( force );
#endif

#if 0
            //printf( "Got force %f\n", force );
            force -= omega2 * displacement;
            force -= damping * velocity;

            accelOld = accel;
            accel = force;

            velocity += ( H / 2.0 ) * accel;
            velocity += ( H / 2.0 ) * accelOld;

            displacement += H * velocity;
            displacement += ( H2 / 2.0 ) * accel;
#endif
#if 0
            velocity += force * H;
            velocity -= omega2 * H * displacement;
            velocity /= vCoefficient;

            displacement += H * velocity;
#endif

            // Try Newmark-beta with beta = 0
            displacement += H * velocity;
            displacement += H2 * accel / 2.0;

            velocity += H * accel / 2.0;
            velocity += ( H / 2.0 ) * ( force - omega2 * displacement );
            velocity = velocity / newmarkVCoefficient;

            accel = force - damping * velocity - omega2 * displacement;

#if 0
            if ( abs(velocity) > 0.0 )
            {
                cout << "Non zero velocity!!!!" << endl;
            }
#endif

            //printf( "Adding sound\n" );
            if ( useDelay ) {
                if ( !buf.add_sound_sample(
                            // Sample value
                            mp_transfer->transfer_norm( ticktime ) * displacement,
                            // Sample time
                            ticktime + delay ) )
                {
                    printf( "Breaking at time %f\n", ticktime );
                    break;
                }
            }
            else {
                if ( !buf.add_sound_sample(
                            mp_transfer->transfer_norm( ticktime ) * displacement ) )
                    //mp_transfer->transfer_norm( ticktime ) * velocity / f ) )
                {
                    printf( "Breaking at time %f\n", ticktime );
                    break;
                }
            }
            //printf( "Done\n" );
        }

#ifdef DEBUG
        char forceProfileFilename[ 1024 ];
        sprintf( forceProfileFilename, "debug_force_profile_%d.vector", mid );

        writeVector( forceProfileFilename, forceProfile );
#endif

#if 0
        if ( omega[mid] < 125.6637 || omega[mid] > CUTTING_OMEGA ) continue;
        printf("freq[%d] = %lf\n", mid, omega[mid]*0.5*M_1_PI);

        // initialize transfer for this mode
        mp_transfer->init_interp(mid);

        const double TSBEGIN = fmax( 0., mp_hertzImpulseSeries->startTime() );
        const double TSEND   = mp_hertzImpulseSeries->endTime();

        mp_hertzImpulseSeries->initInterp( mid );

        double lastm1 = 0, lastm2 = 0;
        //// IIR filter
        double theta = omegad[mid] * H;
        double eps   = exp(-xi[mid] * omega[mid] * H);
        double gamma = asin(xi[mid]);
        double c1 = 2 * eps * cos(theta),
               c2 = - M_SQR(eps),
               c3 = 2*(eps*cos(theta+gamma) + c2*cos(2*theta+gamma))
                   / (3*omega[mid]*omegad[mid]);

        buf.begin_sound_gen(TSBEGIN);
        double ticktime;
        for(int tick = 0;(ticktime = TSBEGIN + tick*H) < TSEND;++ tick)
        {
            //double u = c1*lastm1 + c2*lastm2 + c3*mp_imp->get_impulse(ticktime);
            double u = c1*lastm1 + c2*lastm2
                + c3*mp_hertzImpulseSeries->getForce(ticktime);
            if ( !buf.add_sound_sample(mp_transfer->transfer_norm(ticktime)*u) )
                break;
            lastm2 = lastm1;
            lastm1 = u;
        }
#endif
    }
    return SUCC_RETURN;
}

/*
 * Generate stereo sound
 */
template<class TModalAnalysis, 
    class TImpulseEval, 
    class TTransferEval>
    int RigidSoundObj<TModalAnalysis, TImpulseEval,
    TTransferEval>::generate_stereo_snd(TStereoBuf& buf)
{
    using namespace std;

    if ( !mp_imp->has_impulses() ) return SUCC_RETURN;

    const double H = 1. / (double)buf.sample_rate();
    const int NUM_MODES = mp_modal->num_modes();
    const vector<double>& omega  = mp_modal->omega();
    const vector<double>& omegad = mp_modal->omegad();
    const vector<double>& xi     = mp_modal->xi();
    const vector<double>& gwide  = mp_modal->gaussian_wide();

    for(int mid = 0;mid < NUM_MODES;++ mid)
    {
        if ( omega[mid] < 125.6637 || omega[mid] > CUTTING_OMEGA ) continue;
        printf("freq[%d] = %lf\n", mid, omega[mid]*0.5*M_1_PI);

        // initialize transfer for this mode
        mp_transfer->init_interp(mid);

        const double WIDE    = fmin(gwide[mid], 25.*H);
        const double TSBEGIN = fmax(0., mp_imp->begin_time() - 4*WIDE);
        const double TSEND   = mp_imp->end_time() - log(1E-5)/(xi[mid] * omega[mid]);
        mp_imp->init_interp(mid, WIDE);

        double lastm1 = 0, lastm2 = 0;
        //// IIR filter
        double theta = omegad[mid] * H;
        double eps   = exp(-xi[mid] * omega[mid] * H);
        double gamma = asin(xi[mid]);
        double c1 = 2 * eps * cos(theta),
               c2 = - M_SQR(eps),
               c3 = 2*(eps*cos(theta+gamma) + c2*cos(2*theta+gamma)) / (3*omega[mid]*omegad[mid]);

        buf.begin_sound_gen(TSBEGIN);
        double ticktime, ltf, rtf;
        for(int tick = 0;(ticktime = TSBEGIN + tick*H) < TSEND;++ tick)
        {
            double u = c1*lastm1 + c2*lastm2 + c3*mp_imp->get_impulse(ticktime);
            mp_transfer->stereo_transfer_norm(ticktime, ltf, rtf);
            if ( !buf.add_sound_sample(ltf*u, rtf*u) ) break;
            lastm2 = lastm1;
            lastm1 = u;
        }
    }
    return SUCC_RETURN;
}

/////////////////////////////////////////////////////////////////////
/*
 * Load the impulse recording from file, and add impulses to the corresponding
 * rigid sound objects.
 * The format of the impulse recording file is
 * <timestamp>  <object id(0-based)>  <vertex ID>  <impulse.x>  <impulse.y>
 * <impulse.z>  <T/S>
 * 
 * The last letter "T/S" indicates what kind of vertex ID is used. "T" means
 * the vertex ID is 
 * the id in tetrahedron mesh, whereas "S" means the vertex ID is the id in
 * surface triangle mesh.
 */
template<class TModalAnalysis, 
    class TImpulseEval, 
    class TTransferEval>
static void load_rigid_snd_impulses(
        const char* file, 
        RigidSoundObj<TModalAnalysis, TImpulseEval, TTransferEval>** objlist, 
        int nobj)
{
    typedef RigidSoundObj<TModalAnalysis, TImpulseEval,
            TTransferEval> TRigidSoundObj;
    // get max object id
    int idmax = -1;
    for(int i = 0;i < nobj;++ i) idmax = std::max(idmax, objlist[i]->obj_id());

    // fill the mapping table
    TRigidSoundObj** objs = new TRigidSoundObj*[idmax+1];
    memset(objs, NULL, (idmax+1)*sizeof(TRigidSoundObj*));
    for(int i = 0;i < nobj;++ i) objs[objlist[i]->obj_id()] = objlist[i];

    std::ifstream fin(file);
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot open impulse file: %s\n", file);
        exit(3);
    }

    for(size_t i = 0;i < nobj;++ i) objlist[i]->begin_read_impulse();

    double ts;
    Vector3d fext;
    int objid, vtxid;
    char onsurf;
    char impulseType;
    printf("Read impulse file[%s] ... \n", file);

    int LN = 0;
    //// vtxid is zero-based
    while ( fin >> ts >> objid >> vtxid >> fext.x >> fext.y >> fext.z
            >> onsurf >> impulseType )
    {
        ++ LN;
        if ( objid > idmax || !objs[objid] )
        {
#if 0
            LOGGING_WARNING("objid is out of range when reading "
                    "impulse file [%s] id=%d",
                    file, objid);
#endif
            continue;
        }
        printf("\rload line: %d", LN); fflush(stdout);
        objs[objid]->add_impulse(ts, vtxid, fext, onsurf);
    }

    for(size_t i = 0;i < nobj;++ i) objlist[i]->end_read_impulse();
    delete []objs;
    printf("\nRead impulse file[%s] ... Finished\n", file);
}

/////////////////////////////////////////////////////////////////////
// This version of the function loads rigid impulses with associated
// Hertz time scales
//
/////////////////////////////////////////////////////////////////////
template<class TModalAnalysis, 
    class TImpulseEval, 
    class TTransferEval>
static void load_rigid_snd_impulses(
        const char* file, 
        const char *impulseScalePrefix,
        RigidSoundObj<TModalAnalysis, TImpulseEval, TTransferEval>** objlist, 
        int nobj)
{
    typedef RigidSoundObj<TModalAnalysis, TImpulseEval,
            TTransferEval> TRigidSoundObj;

    char                     buf[ 1024 ];

    // get max object id
    int idmax = -1;
    for(int i = 0;i < nobj;++ i) idmax = std::max(idmax, objlist[i]->obj_id());

    // fill the mapping table
    TRigidSoundObj** objs = new TRigidSoundObj*[idmax+1];
    memset(objs, NULL, (idmax+1)*sizeof(TRigidSoundObj*));
    for(int i = 0;i < nobj;++ i) objs[objlist[i]->obj_id()] = objlist[i];

    bool                     loadedBinaryImpulses = true;
    char                     impFileName[ 1024 ];

    for ( int i = 0; i < nobj; i++ ) {
        sprintf( impFileName, "%s_%s_%d.dat", impulseScalePrefix, file,
                objlist[i]->obj_id() );

        if ( !objlist[ i ]->load_binary_impulses( impFileName ) ) {
            loadedBinaryImpulses = false;
            break;
        }
    }

    if ( !loadedBinaryImpulses )
    {
        // Clear any previously loaded impulses, since we have to start
        // from scratch
        for ( int i = 0; i < nobj; i++ ) {
            objlist[ i ]->hertz_series()->clearImpulses();
        }

        std::ifstream fin(file);
        if ( fin.fail() )
        {
            PRINT_ERROR("Cannot open impulse file: %s\n", file);
            exit(3);
        }

        // Open the impulse scale files
        sprintf( buf, "%s_plane.txt", impulseScalePrefix );
        std::ifstream planeScalesIn( buf );
        if ( planeScalesIn.fail() )
        {
            PRINT_ERROR( "Cannot open scale file %s\n", buf );
            exit(3);
        }

        sprintf( buf, "%s_object.txt", impulseScalePrefix );
        std::ifstream objectScalesIn( buf );
        if ( objectScalesIn.fail() )
        {
            PRINT_ERROR( "Cannot open scale file %s\n", buf );
            exit(3);
        }

#if 0
        for(size_t i = 0;i < nobj;++ i) objlist[i]->begin_read_impulse();
#endif

        double ts;
        Vector3d fext;
        int objid, vtxid;
        char onsurf;
        char impulseType;
        printf("Read impulse file[%s] ... \n", file);

        double impulseLength;
        double impulseScale;

        int planeLine = 0;
        int objLine = 0;

        printf( "Loading scaled impulses\n" );
        printf( "idmax = %d\n", idmax );

        int LN = 0;
        //// vtxid is zero-based
        while ( fin >> ts >> objid >> vtxid >> fext.x >> fext.y >> fext.z
                >> onsurf >> impulseType )
        {
            ++ LN;
#if 0
            if ( objid > idmax || !objs[objid] )
            {
#if 0
                LOGGING_WARNING("objid is out of range when reading "
                        "impulse file [%s] id=%d",
                        file, objid);
#endif
                continue;
            }
#endif
            //printf("\rload line: %d", LN); fflush(stdout);

            if ( impulseType == 'P' )
            {
                // This impulse should be paired with another object.
                // Read the time scale for this pairwise impulse
                if ( objectScalesIn.eof() )
                {
                    printf( "Error reading object scales file\n" );
                }
                objectScalesIn >> impulseLength >> impulseScale;
                objLine += 1;

                int                tetID;

                if ( objid <= idmax && objs[ objid ] )
                {
                    tetID = objs[ objid ]->impulse_interp()->surfid( vtxid, onsurf );

                    if ( impulseLength > 0.0 )
                    {
                        objs[ objid ]->add_impulse( ts, tetID, fext,
                                impulseLength, impulseScale );
                    }
                }

                // This is a pair of impulses, so read the next one as well
                if ( fin.eof() )
                {
                    PRINT_ERROR( "End of impulse file" );
                    exit(3);
                }

                fin >> ts >> objid >> vtxid >> fext.x >> fext.y >> fext.z
                    >> onsurf >> impulseType;

                if ( objid <= idmax && objs[ objid ] )
                {
                    tetID = objs[ objid ]->impulse_interp()->surfid( vtxid, onsurf );

                    if ( impulseLength > 0.0 )
                    {
                        objs[ objid ]->add_impulse( ts, tetID, fext,
                                impulseLength, impulseScale );
                    }
                }
            }
            else
            {
                // This is a single impulse with a constraint (eg. a ground
                // plane
                if ( planeScalesIn.eof() )
                {
                    printf( "Error reading plane scales file\n" );
                }
                planeScalesIn >> impulseLength >> impulseScale;
                planeLine += 1;

                int                tetID;

                if ( objid <= idmax && objs[ objid ] )
                {
                    tetID = objs[ objid ]->impulse_interp()->surfid( vtxid, onsurf );

                    if ( impulseLength > 0.0 )
                    {
                        objs[ objid ]->add_impulse( ts, tetID, fext,
                                impulseLength, impulseScale );
                    }
                }
            }
        }

        if ( !objectScalesIn.eof() )
        {
            printf( "Big problem with the object scale file: %d\n", objLine );

            char tmpBuf[ 4096 ];
            objectScalesIn.getline( tmpBuf, 4096 );
            printf( "Next line: %s\n", tmpBuf );
        }

        if ( !planeScalesIn.eof() )
        {
            printf( "Big problem with the plane scale file: %d\n", planeLine );
            char tmpBuf[ 4096 ];
            planeScalesIn.getline( tmpBuf, 4096 );
            printf( "Next line: %s\n", tmpBuf );
        }

        // Save loaded impulses
        for ( int i = 0; i < nobj; i++ ) {
            sprintf( impFileName, "%s_%s_%d.dat", impulseScalePrefix, file,
                    objlist[i]->obj_id() );

            objlist[ i ]->hertz_series()->save( impFileName );
        }
    }

#if 0
    for(size_t i = 0;i < nobj;++ i) objlist[i]->end_read_impulse();
#endif

    for ( int i = 0; i < nobj; i++ )
    {
        objlist[ i ]->hertz_series()->init();
    }

    delete []objs;
    printf("\nRead impulse file[%s] ... Finished\n", file);
}

#endif

