#include <string>
#include <memory>
#include <Eigen/Dense>
#include <qapplication.h>
#include <wavesolver/SimWorld.h>
#include <ui/FDTD_AcousticSimulator_Viewer.h>
#include <ui/FDTD_AcousticSimulator_Widget.h>
#include <boost/program_options.hpp>

const std::string STR_NOT_EXISTS = "NOT_EXIST";
//##############################################################################
//##############################################################################
int GUI_Run(int argc, char **argv,
            const std::string &xmlFile, const uint &indTimeChunks,
            const uint &previewSpeed);
int TUI_Run(int argc, char **argv,
            const std::string &xmlFile, const uint &indTimeChunks);
void Run_Chunks_Analysis(const std::string &xmlFile);

//##############################################################################
//##############################################################################
int main(int argc, char** argv)
{
    // enable multi-threading with Eigen: https://eigen.tuxfamily.org/dox/TopicMultiThreading.html
    Eigen::initParallel();
    //
    std::string xmlFile, cnkFile;
    uint preview_speed;
    bool nogui;
    uint index = 0u;
    //
    namespace po = boost::program_options;
    try
    {
        po::options_description opt("Options");
        opt.add_options()("help,h", "display help information");
        opt.add_options()
            ("config,c"             , po::value<std::string>(&xmlFile)->required()
                                    , "Configuration xml file")
            ("run_chunks_analysis,a", po::bool_switch()->default_value(false)
                                    , "Analyze time-parallelization chunk sizes")
            ("nogui,n"              , po::bool_switch()->default_value(false)
                                    , "Don't run GUI (if no GUI, number_of_timesteps is needed in settings)")
            ("preview,p"            , po::value<uint>(&preview_speed)->default_value(0)
                                    , "preview rigid body motion")
            ("index,i"              , po::value<uint>(&index)->default_value(0)
                                    , "Time chunk index for time parallelism");
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, opt), vm);
        po::notify(vm);
        if (vm.count("help"))
        {
            std::cout << opt << "\n";
            return 1;
        }
        // assign from parsed commands
        nogui = vm["nogui"].as<bool>();
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

    // if run chunks analysis then skip solver
    if (cnkFile != STR_NOT_EXISTS)
    {
        Run_Chunks_Analysis(xmlFile);
        return 0;
    }

    if (nogui)
        return TUI_Run(argc, argv, xmlFile, index);
    return GUI_Run(argc, argv, xmlFile, index, preview_speed);
}

//##############################################################################
//##############################################################################
int GUI_Run(int argc, char **argv, const std::string &xmlFile,
            const uint &indTimeChunks, const uint &previewSpeed)
{
    QApplication application(argc,argv);

    ImpulseResponseParser_Ptr parser = std::make_shared<ImpulseResponseParser>(xmlFile);
    SimWorld_UPtr world(new SimWorld());
    world->Build(parser, indTimeChunks);
    std::shared_ptr<FDTD_AcousticSimulator_Viewer> viewer =
        std::make_shared<FDTD_AcousticSimulator_Viewer>(std::move(world));
    if (previewSpeed > 0)
        viewer->SetPreviewSpeed(previewSpeed);
    FDTD_AcousticSimulator_Widget widget(viewer);
    widget.setWindowTitle("FDTD Acoustic Simulator Widget");
    widget.show();
    viewer->setWindowTitle("FDTD Acoustic Simulator Viewer");
    viewer->show();

    return application.exec();
}

//##############################################################################
//##############################################################################
int TUI_Run(int argc, char **argv, const std::string &xmlFile, const uint &indTimeChunks)
{
    ImpulseResponseParser_Ptr parser = std::make_shared<ImpulseResponseParser>(xmlFile);
    SimWorld_UPtr world(new SimWorld());
    world->Build(parser, indTimeChunks);
    const int numberSteps = world->GetSolverSettings()->numberTimeSteps;
    for (int ii=0; ii<numberSteps; ++ii)
    {
        world->StepWorld();
    }
    return (numberSteps >=0 ? 0 : -1);
}

//##############################################################################
//##############################################################################
void Run_Chunks_Analysis(const std::string &xmlFile)
{
    std::cout << "\n[Run_Chunks_Analysis]\n";
    std::cout << "==========================================================\n";
    std::cout << "==========================================================\n";
    ImpulseResponseParser_Ptr parser =
        std::make_shared<ImpulseResponseParser>(xmlFile);
    SimWorld_UPtr world(new SimWorld());
    world->Build(parser);
    auto param = parser->GetChunkPartitionParam();
    world->RunChunksAnalysis(param);
    std::cout << "==========================================================\n";
    std::cout << "==========================================================\n";
}
