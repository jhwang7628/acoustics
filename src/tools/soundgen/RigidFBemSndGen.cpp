#include <iostream>
#include <vector>
#include <string>
#include <libconfig.h++>

#include "utils/macros.h"
#include "utils/print_msg.h"
#include "sndgen/FBemTransferInterp.h"
//#include "tests/FixedTransferInterp.h"          // $TEST$
#include "sndgen/RigidModal.h"
#include "sndgen/ImpulseInterp.hpp"
#include "sndgen/RigidSoundObj.hpp"
#include "sndgen/SoundWriter.hpp"

using namespace std;

typedef ImpulseInterp<RigidModal>           TImpulseInterp;
typedef RigidSoundObj<RigidModal, 
                      TImpulseInterp,
                      FBemTransferInterp>   TRigidSoundObj;
//typedef RigidSoundObj<RigidModal,               // $TEST$
//                      TImpulseInterp,
//                      FixedTransferInterp>   TRigidSoundObj;

static double TOT_TIME = 5;
static int SND_RATE = 44100;                      // generated sound rate (Hz)

static vector<TRigidSoundObj*>  sndObjs;
static SoundBuffer              sndBuf;
static SoundBuffer              tsSndBuf;
//static StereoSndBuffer          sndBuf;       // $STEREO_SND$

static map<int, Vector3d>       objectMassCenters;

static void clean()
{
    for(size_t i = 0;i < sndObjs.size();++ i) 
    {
        delete sndObjs[i]->transfer_interp();
        delete sndObjs[i]->impulse_interp();
        delete sndObjs[i]->modal_analysis();
        delete sndObjs[i];
    }
}

static bool load_mass_centers(const char* file)
{
  ifstream                   mcFile( file );

  int                        objID;
  double                     mX, mY, mZ;

  if ( !mcFile.good() ) {
    cerr << "Error opening " << file << endl;
    return false;
  }

  while ( !mcFile.eof() ) {
    mcFile >> objID >> mX >> mY >> mZ;

    printf( "Setting object %d mass center to (%f, %f, %f)\n",
            objID, mX, mY, mZ );

    objectMassCenters[ objID ] = Vector3d( mX, mY, mZ );
  }
}

static int load_config(const char* file, bool &foundImpulseScale)
{
    using namespace libconfig;

    printf("Load configure file: %s\n", file);
    Config cfg;
    try
    {
        string     impulseScalePrefix;

        foundImpulseScale = false;

        cfg.readFile(file);
        double simstep;
        if ( !cfg.lookupValue("sndrate", SND_RATE) ||
             !cfg.lookupValue("sim_step", simstep) ||
             // timestep size used in simulation
             !cfg.lookupValue("sound_len", TOT_TIME) )
        {
            PRINT_ERROR("Cannot find 'sndrate' or 'sim_step' or 'sound_len' "
                        "in configure file: %s\n", file);
            return ERROR_RETURN;
        }

        string massCenterFile;
        if ( !cfg.lookupValue("mass_center_file", massCenterFile) ) {
          cerr << "No mass centers found" << endl;
        } else {
          load_mass_centers( massCenterFile.c_str() );
        }

        Setting& ss1 = cfg.lookup("sndobj");
        int lenss1 = ss1.getLength();
        printf("%d Sound objs will be loaded\n", lenss1);

        if ( cfg.lookupValue("scale_prefix", impulseScalePrefix) )
        {
          foundImpulseScale = true;
        }

        sndObjs.resize(lenss1);
        for(int i = 0;i < lenss1;++ i)
        {
            printf( "Loading object %d\n", i );
            RigidModal* pmodel = new RigidModal(ss1[i]);
            printf( "Built rigid modal\n" );

            Vector3d *massCenter = NULL;

            if ( objectMassCenters.size() > 0 ) {
              int objID = pmodel->id();

              // We should be able to find an entry for this object
              if ( objectMassCenters.find( objID ) == objectMassCenters.end() )
              {
                cerr << "Error; found no mass center for object "
                     << objID << endl;
                exit(3);
              }

              massCenter = &( objectMassCenters[ objID ] );
            }

            TImpulseInterp* pimp = new TImpulseInterp(ss1[i], simstep);
            printf( "Built impulseInterp\n" );
            FBemTransferInterp* ptrans
              = new FBemTransferInterp(ss1[i], pmodel->num_modes(),
                                       pmodel->omega(), pmodel->id(),
                                       massCenter);
            printf( "Built transfer interp\n" );
#if 0
            FixedTransferInterp* ptrans
              = new FixedTransferInterp(ss1[i], pmodel->num_modes(),
                                         pmodel->omega());    // $TEST$
#endif
            if ( foundImpulseScale )
            {
              HertzImpulseSeries<RigidModal>  *hertzImpulseSeries;

              hertzImpulseSeries = new HertzImpulseSeries<RigidModal>( pmodel );
              printf( "Built hertz time series\n" );

              sndObjs[i] = new TRigidSoundObj(ss1[i], pmodel, pimp,
                                              hertzImpulseSeries, ptrans, 0);

              cout << "Generating hertz sound object" << endl;
            }
            else
            {
              sndObjs[i] = new TRigidSoundObj(ss1[i], pmodel, pimp, ptrans, 0);
              //(double)rand()*0.1/(double)RAND_MAX);
            }
        }

        string impfile;
        if ( !cfg.lookupValue("impulse_file", impfile) )
        {
            PRINT_ERROR("Cannot find 'impulse_file' in configure file: %s\n",
                        file);
            return ERROR_RETURN;
        }

        if ( foundImpulseScale )
        {
          load_rigid_snd_impulses(impfile.c_str(),
                                  impulseScalePrefix.c_str(),
                                  (TRigidSoundObj**)&sndObjs[0],
                                  lenss1);
        }
        else
        {
          load_rigid_snd_impulses(impfile.c_str(),
                                  (TRigidSoundObj**)&sndObjs[0],
                                  lenss1);
        }

        return SUCC_RETURN;
    }
    catch (const SettingException& e)
    {
        fprintf(stderr, 
                "Error occured when reading configure file at %s: %s\n",
                e.getPath(), e.what());
        return ERROR_RETURN;
    }
    catch (ParseException& e)
    {
        fprintf(stderr, "ParseException: %s %s at Line %d\n",
                e.getError(), e.what(), e.getLine());
        return ERROR_RETURN;
    }
}

int main(int argc, char* argv[])
{
#if 0
    if ( argc != 3 )
    {
        PRINT_ERROR("Invalid arguments!\n");
        printf("Usage: %s <config_file> <output.wav>\n", argv[0]);
        return 1;
    }
#endif
    if ( argc != 2 )
    {
        PRINT_ERROR("Invalid arguments!\n");
        printf("Usage: %s <config_file>\n", argv[0]);
        return 1;
    }

    bool                     foundImpulseScale = false;
    bool                     useDelay;

    if ( _FAILED(load_config(argv[1], foundImpulseScale)) )
    {
        PRINT_ERROR("Fail to load the configure file: %s\n", argv[1]);
        return 1;
    }

    useDelay = ( objectMassCenters.size() > 0 );

    sndBuf.init(SND_RATE, TOT_TIME);
    tsSndBuf.init(SND_RATE, TOT_TIME);
    for(size_t i = 0;i < sndObjs.size();++ i)
    {
        printf("Generate sound for obj #%d\n", sndObjs[i]->obj_id());
        if ( foundImpulseScale )
        {
          printf( "Generating hertz sound\n" );
          sndObjs[i]->timestep_hertz_sound(tsSndBuf, useDelay);
          //abort();
#if 0
          sndObjs[i]->generate_hertz_sound(sndBuf);
#endif
        }
        else
        {
          sndObjs[i]->generate_sound(sndBuf);
        }
        //sndObjs[i]->generate_stereo_snd(sndBuf);       // $STEREO_SND$
    }
    printf( "Done generating object sounds\n" );

#if 0
    sndBuf.write_to_raw( "test.vector" );
#endif
    if ( objectMassCenters.size() > 0 ) {
      tsSndBuf.write_to_raw( "test_delay.vector" );
    }
    else {
      tsSndBuf.write_to_raw( "test.vector" );
    }

#if 0
    if ( _FAILED(sndBuf.normalize()) )
        PRINT_WARNING("Fail to normalize the sound buffer\n");
#endif

#if 0
    if ( _FAILED(tsSndBuf.normalize() ) )
        PRINT_WARNING("Fail to normalize the time step sound buffer\n");
#endif

#if 0
    char tsFile[ 1024 ];
    sprintf( tsFile, "%s", argv[ 2 ] );

    //// write to wav file
    printf("Writing sound into file: RATE=%d", SND_RATE);
#endif
#if 0
    SingleChannelWavWriter writer(argv[2], SND_RATE, 32);
#endif
#if 0
    printf( "Writing to %s\n", tsFile );
    SingleChannelWavWriter timeStepWriter(tsFile, SND_RATE, 32);
    //StereoWavWriter writer(argv[2], SND_RATE, 32);    // $STEREO_SND$
#endif

#if 0
    writer.set_normalized<double>(true);
    writer.write<double>(sndBuf.data(), sndBuf.num_frames());
#endif

#if 0
    printf( "Normalizing\n" );
    timeStepWriter.set_normalized<double>(true);
    timeStepWriter.write<double>(tsSndBuf.data(), tsSndBuf.num_frames());
#endif

    clean();
    return 0;
}

