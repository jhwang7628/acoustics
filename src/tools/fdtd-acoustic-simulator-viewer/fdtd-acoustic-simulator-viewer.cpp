#include <string>
#include <memory>
#include <qapplication.h>
#include <ui/FDTD_AcousticSimulator_Viewer.h>

int main(int argc, char** argv)
{
    std::string xmlFile("/home/jui-hsien/code/acoustics/src/tools/unit_testing/test_FDTD_RigidObject.xml"); 
    if (argc == 2) 
        xmlFile = std::string(argv[1]);

    QApplication application(argc,argv);

    FDTD_AcousticSimulator_Viewer viewer(xmlFile);
    viewer.setWindowTitle("FDTD Acoustic Simulator Viewer");
    viewer.show();

    return application.exec();
}
