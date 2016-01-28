#ifndef REPLAY_MULTIOBJ_VIEWER
#   define REPLAY_MULTIOBJ_VIEWER

#include <QGLViewer/qglviewer.h>
#include <map>

#include <geometry/TriangleMesh.hpp>
#include <linearalgebra/Quaternion.hpp>
#include <rigid/LSCollisionRigidBody.hpp>
#include <visual/RigidBodyRenderer.hpp>
#include <io/TetMeshReader.hpp>

class MultiViewer : public QGLViewer
{
    Q_OBJECT

    public slots:
        void snapshot(bool);

    public:
        typedef FixVtxTetMesh<double>                   TTetMesh;
        typedef LSCollisionRigidBody<double, TTetMesh>  TRigidBody;

    public:
        MultiViewer(const std::string& disp, const std::string& file):
            displaceFile_(disp), configFile_(file), inUse_(NULL),
            autoSnapshot_(false), showWire_(false)
        {
          load_config();
        }

        ~MultiViewer()
        { 
            delete []inUse_; 
        }

    protected:
        void draw();
        void init();
        void animate();
        void keyPressEvent(QKeyEvent* e);

    private:
        void load_config();
#if 0
        void load_mesh(int id);
#endif

    private:
        struct Pos
        {
            int                 id;
            Vector3<double>     translate;
            Quaternion<double>  rotate;

            Pos() { }
            Pos(int id, const Vector3<double>& t, const Quaternion<double>& r):
                    id(id), translate(t), rotate(r) { }
        };

        std::string                 displaceFile_;
        std::string                 configFile_;
        int                         tsidx_;
        int                         lastts_;
        int                         numObjs_;
        bool*                       inUse_;
        bool                        autoSnapshot_;
        bool                        showWire_;

        std::map< int, TRigidBody* >                            bodies_;
        std::vector< std::pair< double, std::vector<Pos> > >    pos_;
        RigidBodyRenderer                                       renderer_;
};

#endif

