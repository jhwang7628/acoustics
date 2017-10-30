#include <iostream>
#include <libconfig.h++>
#include "utils/print_msg.h"
#include "io/TetMeshReader.hpp"
#include "rigid/Simulator.h"

#include <boost/bind.hpp>

using namespace std;


#define CONFIG_LOOKUP_HELPER(c,tag,variable) \
    c.lookupValue(tag,variable);

#define CONFIG_LOOKUP_HELPER_DEFAULT(c,tag,variable,defaultValue) \
    if(!c.lookupValue(tag,variable)) variable=defaultValue;

#define CONFIG_LOOKUP_HELPER_ASSERT(c,tag,variable) \
    assert(c.lookupValue(tag,variable))

#define CONFIG_LOOKUP_HELPER_FLAG(c,tag,variable,flag) \
    flag=c.lookupValue(tag,variable);


// ----------------------------------------------------------------------------
struct MatRec
{
    double density;
    double rest;
    double friction;
};

// Animation state function to bind
double sinPulse( double t, double tStart, double tEnd,
                 double period, double amplitude,
                 double offset = 0.0 )
{
  if ( t < tStart )
  {
    return offset;
  }

  if ( t > tEnd )
  {
    t = tEnd;
  }

  double                     signalResult;

  signalResult = t - tStart;
  signalResult /= period;
  signalResult *= 2.0 * M_PI;
  signalResult = sin( signalResult );
  signalResult *= amplitude;

  return signalResult + offset;
}

double cosPulse( double t, double tStart, double tEnd,
                 double period, double amplitude,
                 double offset = 0.0 )
{
  if ( t < tStart || t > tEnd )
    return offset;

  double                     signalResult;

  signalResult = t - tStart;
  signalResult /= period;
  signalResult *= 2.0 * M_PI;
  signalResult = cos( signalResult );
  signalResult *= amplitude;

  return signalResult + offset;
}

double zeroFunc( double t, double offset = 0.0 )
{
  return offset;
}

typedef FixVtxTetMesh<double>           TTetMesh;
typedef RigidBodySimulator::TRigidBody  TRigidBody;

static RigidBodySimulator   sim;
static double               stepSize, timeLen;

static void load_config(const char* file)
{
    using namespace libconfig;

    map<string, MatRec>  mats;
    string txt;
    Config cfg;
    MatRec mat;
    try
    {
        cfg.readFile(file);
        if ( !cfg.lookupValue("simulation.step_size", stepSize) )
        {
            PRINT_ERROR("Cannot find 'step_size' in config file %s\n", file);
            exit(1);
        }
        if ( !cfg.lookupValue("simulation.time_len", timeLen) )
        {
            PRINT_ERROR("Cannot find 'time_len' in config file %s\n", file);
            exit(1);
        }

        Setting& ss = cfg.lookup("materials");
        int sslen = ss.getLength();
        printf("INFO: Read %d materials ...\n", sslen);
        for(int i = 0;i < sslen;++ i)
        {
            if ( !ss[i].lookupValue("name", txt) ||
                 !ss[i].lookupValue("density", mat.density) ||
                 !ss[i].lookupValue("rest_coeff", mat.rest) ||
                 !ss[i].lookupValue("friction_coeff", mat.friction) )
            {
                PRINT_ERROR("Matieral section is incomplete\n");
                exit(1);
            }
            mats[txt] = mat;
        }

        Setting& ssO = cfg.lookup("objs");
        sslen = ssO.getLength();
        TTetMesh* pmesh = NULL;

        map<string, TTetMesh*>  meshtbl;
        vector<TRigidBody*>     rbodies(sslen);
        printf("INFO: Read %d objects ...\n", sslen);
        for(int i = 0;i < sslen;++ i)
        {
            if ( !ssO[i].lookupValue("model", txt) )
            {
                PRINT_ERROR("No model name specified\n");
                exit(1);
            }
            if ( !meshtbl.count(txt) )
            {
                pmesh = new TTetMesh;
                FV_TetMeshLoader_Double::load_mesh(txt.c_str(), *pmesh);
                pmesh->init();
                pmesh->update_surface();
                meshtbl[txt] = pmesh;
            }
            else 
                pmesh = meshtbl[txt];

            if ( !ssO[i].lookupValue("material", txt) )
            {
                PRINT_ERROR("No material specified\n");
                exit(1);
            }
            if ( !mats.count(txt) )
            {
                PRINT_ERROR("Cannot find desired material: %s\n", txt.c_str());
                exit(1);
            }
            const MatRec& matref = mats[txt];
            rbodies[i] = new TRigidBody(i, pmesh, matref.density, 
                                        matref.rest, matref.friction);
            double dx, dy, dz, rot = 0.;
            Quaternion<double> qrot;
            /*
             * Initialize object state
             */
            if ( !ssO[i].lookupValue("dx", dx) ) dx = 0;
            if ( !ssO[i].lookupValue("dy", dy) ) dy = 0;
            if ( !ssO[i].lookupValue("dz", dz) ) dz = 0;
            std::vector<Quat4d> rotations;
            if (ssO[i].exists("rotations"))
            {
                const libconfig::Setting &r = ssO[i]["rotations"];
                double angle, rx, ry, rz;
                for(int j=0; j<r.getLength(); j++)
                {
                    CONFIG_LOOKUP_HELPER_ASSERT(r[j],"angle",angle);
                    CONFIG_LOOKUP_HELPER_ASSERT(r[j],"rx"   ,rx);
                    CONFIG_LOOKUP_HELPER_ASSERT(r[j],"ry"   ,ry);
                    CONFIG_LOOKUP_HELPER_ASSERT(r[j],"rz"   ,rz);
                    rotations.push_back(Quat4d::fromAxisRotD(Vector3d(rx,ry,rz), angle));
                }
            }
            if (rotations.size()>0)
            {
                Quat4d quat = rotations[0];
                for(int j=1;j<rotations.size();j++)
                    quat *= rotations[j];
                rbodies[i]->init_com_rotation(quat.w,quat.v.x,quat.v.y,quat.v.z);
            }
            rbodies[i]->translate(dx, dy, dz);

#if 0
            // FIXME
            qrot = Quat4d::fromAxisRotD(Vector3d(1,0,0), 90.0);
            rbodies[i]->rotate( qrot );
#endif
#if 0
            if ( ssO[i].lookupValue("rx", rot) )
                qrot = Quat4d::fromAxisRotD(Vector3d(1,0,0), rot);
            else if ( ssO[i].lookupValue("ry", rot) )
                qrot = Quat4d::fromAxisRotD(Vector3d(0,1,0), rot);
            else if ( ssO[i].lookupValue("rz", rot) )
                qrot = Quat4d::fromAxisRotD(Vector3d(0,0,1), rot);
            else 
                rot = 0.;
            cerr << "rot = " << rot << std::endl;

            if ( fabs(rot) > 1E-9 )
            {
                Point3d ctr = rbodies[i]->mass_center();
                rbodies[i]->set_position(Point3d(ctr.x+dx, ctr.y+dy, ctr.z+dz), qrot);
            }
            else
                rbodies[i]->translate(dx, dy, dz);
#endif

            int setFixed = 0;
            if ( !ssO[i].lookupValue("fixed", setFixed) ) setFixed = 0;

            if ( setFixed )
            {
              printf( "Fixing an object!\n" );
              rbodies[i]->set_fixed( true );
            }
            else
            {
              int             animate = 0;

              if ( !ssO[i].lookupValue("animate", animate) ) animate = 0;

              if ( animate )
              {
                Point3d pos = rbodies[i]->mass_center();

                const Setting &animateDirection = ssO[i]["animate_direction"];
                double         animateAmplitude;
                double         animateStart;
                double         animateEnd;
                double         animatePeriod;

                if ( !ssO[i].lookupValue("animate_amplitude",
                                         animateAmplitude ) )
                {
                  animateAmplitude = 0.0;
                }
                if ( !ssO[i].lookupValue("animate_start",
                                         animateStart ) )
                {
                  animateStart = 0.0;
                }
                if ( !ssO[i].lookupValue("animate_end",
                                         animateEnd ) )
                {
                  animateEnd = 0.0;
                }
                if ( !ssO[i].lookupValue("animate_period",
                                         animatePeriod ) )
                {
                  animatePeriod = 0.0;
                }

                // Check to see if we want to specify an animation state
                //
                // NOTE: I hacked this in here for the key chain example
                RigidState     *s = new RigidState();

                s->positionSignal[ 0 ] = boost::bind(
                              sinPulse, _1, animateStart, animateEnd,
                              animatePeriod,
                              (double)animateDirection[ 0 ] * animateAmplitude,
                              pos[0] );
                s->positionSignal[ 1 ] = boost::bind(
                              sinPulse, _1, animateStart, animateEnd,
                              animatePeriod,
                              (double)animateDirection[ 1 ] * animateAmplitude,
                              pos[1] );
                s->positionSignal[ 2 ] = boost::bind(
                              sinPulse, _1, animateStart, animateEnd,
                              animatePeriod,
                              (double)animateDirection[ 2 ] * animateAmplitude,
                              pos[2] );

                s->velocitySignal[ 0 ] = boost::bind( cosPulse, _1,
                                  animateStart,
                                  animateEnd,
                                  animatePeriod,
                                  2.0 * M_PI * (double)animateDirection[ 0 ]
                                    * animateAmplitude / animatePeriod,
                                  0.0 );
                s->velocitySignal[ 1 ] = boost::bind( cosPulse, _1,
                                  animateStart,
                                  animateEnd,
                                  animatePeriod,
                                  2.0 * M_PI * (double)animateDirection[ 1 ]
                                    * animateAmplitude / animatePeriod,
                                  0.0 );
                s->velocitySignal[ 2 ] = boost::bind( cosPulse, _1,
                                  animateStart,
                                  animateEnd,
                                  animatePeriod,
                                  2.0 * M_PI * (double)animateDirection[ 2 ]
                                    * animateAmplitude / animatePeriod,
                                  0.0 );

#if 0
                s->positionSignal[ 1 ] = boost::bind( zeroFunc, _1, pos[1] );
                s->positionSignal[ 2 ] = boost::bind( zeroFunc, _1, pos[2] );

                s->velocitySignal[ 1 ] = boost::bind( zeroFunc, _1, 0.0 );
                s->velocitySignal[ 2 ] = boost::bind( zeroFunc, _1, 0.0 );
#endif

                rbodies[i]->specify_state( s );
              }
            }

            if ( !ssO[i].lookupValue("vx", dx) ) dx = 0;
            if ( !ssO[i].lookupValue("vy", dy) ) dy = 0;
            if ( !ssO[i].lookupValue("vz", dz) ) dz = 0;
            rbodies[i]->init_velocity(dx, dy, dz);
            if ( !ssO[i].lookupValue("wx", dx) ) dx = 0;
            if ( !ssO[i].lookupValue("wy", dy) ) dy = 0;
            if ( !ssO[i].lookupValue("wz", dz) ) dz = 0;
            rbodies[i]->init_angular_velocity(dx, dy, dz);
            /* add to rigid simulator */
            sim.add_rigid_body(rbodies[i]);
        }

        //// clean unused tet mesh to save memory
        for(size_t i = 0;i < rbodies.size();++ i)
            rbodies[i]->clean_tet_mesh();
        map<string, TTetMesh*>::iterator end = meshtbl.end();
        for(map<string, TTetMesh*>::iterator it = meshtbl.begin(); it != end;++ it)
            delete it->second;

        printf("INFO: Initialize simulator ...\n");
        sim.init();
    }
    catch (const SettingException e)
    {
        fprintf(stderr, "Error occured when reading configure file at %s: %s\n",
                e.getPath(), e.what());
    }
    catch (ParseException& e)
    {
        fprintf(stderr, 
                "ParseException: %s %s at Line %d\n",
                e.getError(), e.what(), e.getLine());
    }
}

int main(int argc, char* argv[])
{
    if ( argc != 2 )
    {
        PRINT_ERROR("Error: Usage: %s <cfg file>\n", argv[0]);
        return 1;
    }

    load_config(argv[1]);
    while ( sim.time() < timeLen )
    {
        printf( "%03.5f\r", sim.time() );
        fflush(stdout);
        //cout << sim.time() << endl;
        sim.advance(stepSize);
    }
    printf( "\n" );

    cout << "Done!" << endl;
    return 0;
}

