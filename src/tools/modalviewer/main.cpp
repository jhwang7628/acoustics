#include <string>
#include <memory>
#include <qapplication.h>
#include <ui/ModalViewer.h>
#include <wavesolver/FDTD_RigidSoundObject.h>

int main(int argc, char** argv)
{
    const std::string meshFileName("/home/jui-hsien/code/acoustics/work/plate/plate.obj");
    const std::string sdfFilePrefix("/home/jui-hsien/code/acoustics/work/plate/plate.obj.1m.dist");
    const int sdfResolution = 100; 
    std::shared_ptr<FDTD_RigidSoundObject> object = std::make_shared<FDTD_RigidSoundObject>(meshFileName, sdfResolution, sdfFilePrefix); 

    QApplication application(argc,argv);

    ModalViewer viewer(object);
    viewer.setWindowTitle("Modal Viewer");
    viewer.show();

    return application.exec();
}
