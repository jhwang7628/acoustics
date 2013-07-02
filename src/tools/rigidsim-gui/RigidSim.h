#ifndef UI_RIGID_SIM_H
#   define UI_RIGID_SIM_H

#include <QApplication>
#include <QLabel>
#if QT_VERSION >= 0x040000
#   include "ui_RigidSim.h"
#else
#   error Only Qt with the version higher than 4.0 is supported right now
#endif

class XDemo;

class RigidSim : public QMainWindow, private Ui_RigidSim
{
    Q_OBJECT

    public slots:
        void start();
        void step();

        void demo_shaking_piggy();
        void demo_drop_objects();
        void demo_drop_objs_with_fixed();

    public:
        RigidSim():demoId_(DEMO_UNKNOWN), pdemo_(NULL)
        {
            setupUi(this);
            statusbar->addWidget(new QLabel("  Press \"H\" for help   ", statusbar));            

            QObject::connect(actionStart, SIGNAL(triggered()), this, SLOT(start()));
            QObject::connect(actionStep,  SIGNAL(triggered()), this, SLOT(step()));

            QObject::connect(actionDropObjects,  SIGNAL(triggered()), this, SLOT(demo_drop_objects()));
            QObject::connect(actionDropObjsWithFixed, SIGNAL(triggered()), this, SLOT(demo_drop_objs_with_fixed()));
            QObject::connect(actionShakingPiggy, SIGNAL(triggered()), this, SLOT(demo_shaking_piggy()));
        }

        ~RigidSim();

    private:
        enum DEMO_ID
        {
            DEMO_UNKNOWN = 0,
            DEMO_SHAKING_PIGGY,
            DEMO_DROP_OBJS,
            DEMO_DROP_OBJS_WITH_FIXED
        };

        DEMO_ID     demoId_;
        XDemo*      pdemo_;
};

#endif
