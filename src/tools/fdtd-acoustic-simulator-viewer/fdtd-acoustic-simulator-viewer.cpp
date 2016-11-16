#include <string>
#include <memory>
#include <qapplication.h>
#include <ui/FDTD_AcousticSimulator_Viewer.h>
#include <ui/FDTD_AcousticSimulator_Widget.h>
#include <boost/program_options.hpp>


int main(int argc, char** argv)
{
    // 
    std::string xmlFile; 
    uint preview_speed; 
    //
    namespace po = boost::program_options; 
    try
    {
        po::options_description opt("Options"); 
        opt.add_options()("help,h", "display help information"); 
        opt.add_options()
            ("config,c", po::value<std::string>(&xmlFile)->required(), "configuration xml file")
            ("preview_speed,p", po::value<uint>(&preview_speed)->default_value(0), "preview rigid body motion"); 
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, opt), vm);
        po::notify(vm);
        if (vm.count("help"))
        {
            std::cout << opt << "\n"; 
            return 1; 
        }
    }
    catch(std::exception &e)
    {
        std::cerr << "**ERROR** " << e.what() << "\n"; 
        return 1; 
    }
    catch(...)
    {
        std::cerr << "**ERROR** Unknown error!" << "\n";
        return false;
    }

    QApplication application(argc,argv);

    std::shared_ptr<FDTD_AcousticSimulator_Viewer> viewer = std::make_shared<FDTD_AcousticSimulator_Viewer>(xmlFile, preview_speed);
    FDTD_AcousticSimulator_Widget widget(viewer); 
    widget.setWindowTitle("FDTD Acoustic Simulator Widget");
    widget.show();
    viewer->setWindowTitle("FDTD Acoustic Simulator Viewer");
    viewer->show();

    return application.exec();
}
