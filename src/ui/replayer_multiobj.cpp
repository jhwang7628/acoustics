#include <iostream>
#include <fstream>
#include <string>
#include <QKeyEvent>
#include <QApplication>

#if defined(__APPLE__) || defined(MACOSX)
#   include <glut.h>
#elif defined(__linux)
#   include <GL/glut.h>
#else
#   error ERROR Unsupported OS
#endif

#include <libconfig.h++>

#include "replayer_multiobj.h"

#include "utils/print_msg.h"

#undef USE_RECORDER

int frameRate = 100;

using namespace std;

void MultiViewer::snapshot(bool at)
{
    //if ( lastts_ == tsidx_ || tsidx_ < 0 || tsidx_ % 100 ) return;
    if ( lastts_ == tsidx_ || tsidx_ < 0 ) return;
    lastts_ = tsidx_;
    saveSnapshot(at);
}

void MultiViewer::draw()
{
    renderer_.set_color(.15f, .15f, 0.15f);

    if ( tsidx_ >= 0 )
    {
        for(size_t i = 0;i < pos_[tsidx_].second.size();++ i)
        {
            const Pos& p = pos_[tsidx_].second[i];
            renderer_.render(bodies_[p.id], p.translate, p.rotate);
        }
    }

    /////////////////////////////////////////////
    const REAL GD_SIZE = 0.01;
    REAL step = GD_SIZE * 15;
    glColor3f(0.7, 0.7, 0.7);

    REAL d = step;
    for(int i = 0;i < 100;++ i, d += step)
    {
        glBegin(GL_LINE_LOOP);
        glVertex3d(-d, 0, -d);
        glVertex3d( d, 0, -d);
        glVertex3d( d, 0,  d);
        glVertex3d(-d, 0,  d);
        glEnd();
    }
#if 0
    drawAxis();
#endif
}

void MultiViewer::keyPressEvent(QKeyEvent* e)
{
    const int key = e->key();
    const Qt::KeyboardModifiers modifiers = e->modifiers();
    
    if ( key == Qt::Key_Space )
        toggleAnimation();
    else if ( key == Qt::Key_P && modifiers == Qt::NoButton )
    {
        autoSnapshot_ = !autoSnapshot_;
        if ( autoSnapshot_ )
        {
            printf("Auto snapshot: ON\n");
            QObject::connect(this, SIGNAL(drawFinished(bool)), SLOT(snapshot(bool)));
        }
        else
        {
            printf("Auto snapshot: OFF\n");
            QObject::disconnect(this, SLOT(snapshot(bool)));
        }
    }
    else if ( key == Qt::Key_N && modifiers == Qt::NoButton )
    {
        animate();
        updateGL();
    }
    else if ( key == Qt::Key_I && modifiers == Qt::NoButton )
    {
        printf("POSITIONS: \n");
        for(int i = 0;i < pos_[tsidx_].second.size();++ i)
        {
            const Pos& p = pos_[tsidx_].second[i];
            printf("  #%d: [%.15lf,%.15lf,%.15lf] [%.15lf;%.15lf,%.15lf,%.15lf]\n",
                    p.id, p.translate.x, p.translate.y, p.translate.z,
                    p.rotate.w, p.rotate.v.x, p.rotate.v.y, p.rotate.v.z);
        }
    }
    else if ( key == Qt::Key_J && modifiers == Qt::NoButton )
    {
        int nf;
        cout << "# of frames to jump over: ";
        cin >> nf;
        tsidx_ = std::max(0, std::min(tsidx_ + nf, (int)pos_.size()-1));
    }
    else if ( key == Qt::Key_W && modifiers == Qt::NoButton )
    {
        showWire_ = !showWire_;
        if ( showWire_ )
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        else
            glPolygonMode(GL_FRONT, GL_FILL);
    }
    else
        QGLViewer::keyPressEvent(e);
}

struct MatRec
{
    double density;
    double rest;
    double friction;
};

void MultiViewer::load_config()
{
    using namespace libconfig;

    map<string, MatRec>  mats;
    string txt;
    Config cfg;
    MatRec mat;
    try
    {
        cfg.readFile(configFile_.c_str());
#if 0
        if ( !cfg.lookupValue("simulation.step_size", stepSize_) )
        {
            PRINT_ERROR("Cannot find 'step_size' in config file %s\n",
                        configFile_.c_str());
            exit(1);
        }
        if ( !cfg.lookupValue("simulation.time_len", timeLen_) )
        {
            PRINT_ERROR("Cannot find 'time_len' in config file %s\n",
                        configFile_.c_str());
            exit(1);
        }
#endif

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
                //TetMeshLoader_Double::load_mesh(txt.c_str(), *pmesh);
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
             *
             * NOTE: assuming the translation and rotation are all related to the origin
             */
#if 0
            if ( !ssO[i].lookupValue("dx", dx) ) dx = 0;
            if ( !ssO[i].lookupValue("dy", dy) ) dy = 0;
            if ( !ssO[i].lookupValue("dz", dz) ) dz = 0;
            rbodies[i]->translate(dx, dy, dz);

            if ( ssO[i].exists("rot") ) 
            {
                const Setting& rr = ssO[i]["rot"];
                rbodies[i]->init_origin_rotation(rr[0], rr[1], rr[2], rr[3]);
            }

#if 0
            // FIXME
            qrot = Quat4d::fromAxisRotD(Vector3d(1,0,0), 90.0);
            rbodies[i]->rotate( qrot );
#endif

            /*
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
            */

            if ( !ssO[i].lookupValue("vx", dx) ) dx = 0;
            if ( !ssO[i].lookupValue("vy", dy) ) dy = 0;
            if ( !ssO[i].lookupValue("vz", dz) ) dz = 0;
            rbodies[i]->init_velocity(dx, dy, dz);
            /* add to rigid simulator */
#endif
        }

        for ( int i = 0; i < rbodies.size(); i++ )
        {
          printf( "Setting body %d\n", i );
          bodies_[ i ] = rbodies[ i ];
        }

#if 0
        //// clean unused tet mesh to save memory
        for(size_t i = 0;i < rbodies.size();++ i)
            rbodies[i]->clean_tet_mesh();
        map<string, TTetMesh*>::iterator end = meshtbl.end();
        for(map<string, TTetMesh*>::iterator it = meshtbl.begin(); it != end;++ it)
            delete it->second;
#endif
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

#if 0
void MultiViewer::load_mesh(int id)
{
    char filename[128];
    sprintf(filename, tetMshPtn_.c_str(), id);
    TTetMesh* pmesh = new TTetMesh;
    printf("INFO: read mesh %s\n", filename);
    //if ( id < 29 )
    //    TetMeshLoader_Double::load_mesh(filename, pmesh, 0.0085, Vector3d(0.,0.,0.));    // $DEMO: Smash-Wall
    //else
    //    FV_TetMeshLoader_Double::load_mesh(filename, *pmesh);    // $DEMO: Smash-Wall
    FV_TetMeshLoader_Double::load_mesh(filename, *pmesh);
    pmesh->init();

    bodies_[id] = new TRigidBody(id, pmesh, 2700, 0.6, 0.6);
}
#endif

void MultiViewer::animate()
{
    //tsidx_ = std::min(tsidx_ + 1, (int)pos_.size()-1);
    tsidx_ = std::min(tsidx_ + frameRate, (int)pos_.size()-1);
    if ( tsidx_ < 0 ) return;

    cerr << "Frame#: " << tsidx_ << endl;
    memset(inUse_, false, sizeof(bool)*numObjs_);
    for(size_t i = 0;i < pos_[tsidx_].second.size();++ i)
    {
        const Pos& p = pos_[tsidx_].second[i];
        inUse_[p.id] = true;
#if 0
        if ( !bodies_.count(p.id) ) load_mesh(p.id);
#endif
        if ( !bodies_.count(p.id) )
        {
          printf( "ERROR: couldn't find body %d\n", p.id );
          abort();
        }
    }

#if 0
    for(size_t i = 0;i < numObjs_;++ i)
        if ( !inUse_[i] && bodies_.count(i) ) 
        {
            TRigidBody* pb = bodies_[i];
            bodies_.erase(i);
            TTetMesh*  msh = pb->mesh();
            delete msh; delete pb;
        }
#endif
}

static GLfloat GLOBAL_AMBIENT[] = { 0.2f, 0.2f, 0.2f, 1.0f };
static GLfloat SPECULAR_COLOR[] = { 0.1f, 0.1f, 0.1f, 1.0 };

void MultiViewer::init()
{
#ifdef __linux
    int dummy = 0;
    glutInit(&dummy, NULL);
#endif
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, GLOBAL_AMBIENT);
    glShadeModel(GL_SMOOTH);

    const GLfloat ambientLight[] = { 0.f, 0.f, 0.f, 1.0f };
    const GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8f, 1.0f };
    const GLfloat specularLight[] = { 0.1f, 0.1f, 0.1f, 1.0f };
    const GLfloat position[] = { -0.5f, 1.0f, 0.4f, 1.0f };

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 1.);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, SPECULAR_COLOR);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    /////
    resize(1024, 576);
    camera()->setZNearCoefficient(0.0001f);
    camera()->setZClippingCoefficient(100.f);

    setAnimationPeriod(0);
    setSnapshotFormat("PNG");
    setSnapshotFileName("frame");

    //// load the displacement
    ifstream fin(displaceFile_.c_str(), ios::binary);
    double ts;
    int    id;
    Point3<REAL>     translate;
    Quaternion<REAL> rotate;
    numObjs_ = 0;

    fin.read((char *)&ts, sizeof(double));
    while ( !fin.fail() )
    {
        pos_.push_back(make_pair(ts, vector<Pos>()));

        fin.read((char *)&id, sizeof(int));
        while ( !fin.fail() && id>= 0 )
        {
            numObjs_ = std::max(id, numObjs_);
            fin.read((char *)&translate, sizeof(Point3<REAL>));
            fin.read((char *)&rotate, sizeof(Quaternion<REAL>));

            if ( fin.fail() )
            {
                cerr << "WARNING: Data is incomplete" << endl;
                break;
            }

            pos_.back().second.push_back(Pos(id, translate, rotate));
            fin.read((char *)&id, sizeof(int));
        }
        fin.read((char *)&ts, sizeof(double));
    }
    printf("INFO: %d timestep loaded\n", (int)pos_.size());

    fin.close();

    tsidx_ = lastts_ = -2;
    ++ numObjs_;
    inUse_ = new bool[numObjs_];
}

///////////////////////////////////////////////////////////////////////////////
static void usage(const char* cmd)
{
    printf("Usage: %s -d <displace file> -f <surface mesh>\n", cmd);
}

static string configFile;
static string displaceFile;

static void parse_cmd(int argc, char* argv[])
{
    int opt;
    while ( (opt = getopt(argc, argv, "hd:f:r:")) != -1 )
    {
        switch (opt)
        {
            case 'h':
                usage(argv[0]);
                exit(0);
            case 'f':
                configFile = optarg;
                break;
            case 'd':
                displaceFile = optarg;
                break;
            case 'r':
                frameRate = atoi( optarg );
                printf( "Frame rate = %d\n", frameRate );
                break;
        }
    }
}

int main(int argc, char** argv)
{
    parse_cmd(argc, argv);

    QApplication application(argc, argv);

    MultiViewer viewer(displaceFile, configFile);

    viewer.setWindowTitle("replayer");
    viewer.show();
    return application.exec();
}
