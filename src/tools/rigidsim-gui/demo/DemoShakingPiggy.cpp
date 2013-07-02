#include "DemoShakingPiggy.h"
#if defined(__APPLE__) || defined(MACOSX)
#   include <glut.h>
#elif defined(__linux)
#   include <GL/glut.h>
#else
#   error ERROR Unsupported OS
#endif

#include <libconfig.h++>
#include <QGLViewer/qglviewer.h>
#include <map>
#include <string>

#include "io/TetMeshReader.hpp"
#include "utils/print_msg.h"

using namespace std;

// ----------------------------------------------------------------------------
struct MatRec
{
    double density;
    double rest;
    double friction;
};

// ----------------------------------------------------------------------------
DemoShakingPiggy::DemoShakingPiggy(const char* file, QGLViewer* canvas):
        canvas_(canvas)
{
    using namespace libconfig;

    map<string, MatRec>  mats;
    string  txt;
    Config  cfg;
    MatRec  mat;
    try
    {
        cfg.readFile(file);

        if ( !cfg.lookupValue("simulation.step_size", stepSize_) )
        {
            PRINT_ERROR("Cannot find 'step_size' in config file %s\n", file);
            exit(1);
        }
        if ( !cfg.lookupValue("simulation.time_len", timeLen_) )
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
                TetMeshLoader_Double::load_mesh(txt.c_str(), pmesh);
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
            double dx, dy, dz;
            /*
             * NOTE: assuming the translation and rotation are all related to the origin
             */
            if ( !ssO[i].lookupValue("dx", dx) ) dx = 0;
            if ( !ssO[i].lookupValue("dy", dy) ) dy = 0;
            if ( !ssO[i].lookupValue("dz", dz) ) dz = 0;

            if ( !i ) 
            {
                Point3d ctr = rbodies[i]->mass_center();
                rbodies[i]->set_position(Point3d(ctr.x+dx, ctr.y+dy, ctr.z+dz), 
                                         Quat4d::fromAxisRotD(Vector3d(0,0,1), 90.));
                piggy_        = rbodies[i];
                piggyInitPos_ = piggy_->mass_center();
                piggyInitRot_ = piggy_->rotation();
            }
            else
                rbodies[i]->translate(dx, dy, dz);

            /* add to rigid simulator */
            rsim_.add_rigid_body(rbodies[i]);
        }

        //// clean unused tet mesh to save memory
        for(size_t i = 0;i < rbodies.size();++ i)
            rbodies[i]->clean_tet_mesh();
        map<string, TTetMesh*>::iterator end = meshtbl.end();
        for(map<string, TTetMesh*>::iterator it = meshtbl.begin(); it != end;++ it)
            delete it->second;

        printf("INFO: Initialize simulator ...\n");
        rsim_.init();
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

int DemoShakingPiggy::start()
{
    const double OMEGA = 4*M_PI;
    const double A     = 0.2;
    
    double tm;
    do {
        rsim_.advance(stepSize_);
        // project the position and vel. of the piggy
        // to make the piggy shake
        tm = rsim_.time();
        piggy_->set_position(
                Point3d(piggyInitPos_.x,
                        piggyInitPos_.y + A*sin(OMEGA*tm),
                        piggyInitPos_.z),
                piggyInitRot_);
        piggy_->set_velocity(
                Vector3d(0., OMEGA*A*cos(OMEGA*tm), 0.),
                Vector3d(0., 0., 0.));
        canvas_->updateGL();
    } while ( tm < timeLen_ );

    printf("Rigid simulation is done\n");
    return 0;
}

int DemoShakingPiggy::step()
{
    PRINT_WARNING("step() is now supported\n");
    return 0;
}

void DemoShakingPiggy::draw()
{
    rbRender_.set_color(0.f, 1.f, 0.4f);
    const std::vector<TRigidBody*>& objs = rsim_.rigid_bodies();
    for(size_t i = 0;i < objs.size();++ i)
    {
        rbRender_.render(objs[i]);
    }
}

