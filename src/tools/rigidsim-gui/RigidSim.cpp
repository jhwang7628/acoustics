#include "RigidSim.h"
#include <QFileDialog>

#include "logging/logging.h"
#include "demo/DemoShakingPiggy.h"
#include "demo/DemoDropObjs.h"
#include "demo/DemoDropObjsWithFixed.h"

RigidSim::~RigidSim()
{   delete pdemo_; }

void RigidSim::start()
{   if ( pdemo_ ) pdemo_->start(); }

void RigidSim::step()
{   if ( pdemo_ ) pdemo_->step(); }

/* 
 * demo: shaking the piggy bank
 */
void RigidSim::demo_shaking_piggy()
{
    if ( demoId_ != DEMO_UNKNOWN ) return;
    QString file = QFileDialog::getOpenFileName(this, "Select the configure file",
            ".", "config file (*.cfg);;All files (*)");
    if ( file.isEmpty() ) return;

    demoId_ = DEMO_SHAKING_PIGGY;
    pdemo_  = new DemoShakingPiggy(file.toStdString().data(), canvas);
    canvas->pdemo_ = pdemo_;
    canvas->updateGL();
}

void RigidSim::demo_drop_objects()
{
    if ( demoId_ != DEMO_UNKNOWN ) return;
    QString file = QFileDialog::getOpenFileName(this, "Select the configure file",
            ".", "config file (*.cfg);;All files (*)");
    if ( file.isEmpty() ) return;

    demoId_ = DEMO_DROP_OBJS;
    pdemo_  = new DemoDropObjs(file.toStdString().data(), canvas);
    canvas->pdemo_ = pdemo_;
    canvas->updateGL();
}

void RigidSim::demo_drop_objs_with_fixed()
{
    if ( demoId_ != DEMO_UNKNOWN ) return;
    QString file = QFileDialog::getOpenFileName(this, "Select the configure file",
            ".", "config file (*.cfg);;All files (*)");
    if ( file.isEmpty() ) return;
    demoId_ = DEMO_DROP_OBJS_WITH_FIXED;
    pdemo_ = new DemoDropObjsWithFixed(file.toStdString().data(), canvas);
    canvas->pdemo_ = pdemo_;
    canvas->updateGL();
}
//=============================================================================

int main(int argc, char* argv[])
{
    LoggingManager::instance().set_logging_level(LOG_INFO);

    QApplication    application(argc, argv);
    RigidSim        mainwnd;
    mainwnd.show();

    return application.exec();
}

