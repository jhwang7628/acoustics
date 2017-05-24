#include "wavesolver/SimWorld.h" 

void TestBuildWorld(const std::string &xmlName)
{
    // build parser
    ImpulseResponseParser_Ptr parser = std::make_shared<ImpulseResponseParser>(xmlName); 

    // build world
    SimWorld world; 
    world.Build(parser); 

}

//##############################################################################
int main(int argc, char ** argv)
{
    if (argc < 2) 
    {
        std::cerr << "**Usage: " << argv[0] << " <config_xml> [steps]\n"; 
        return -1; 
    }
    const std::string xmlName = std::string(argv[1]);
    const int steps = (argc == 3 ? atoi(argv[2]) : -1); 
    {
        boost::timer::auto_cpu_timer timer;
        TestBuildWorld(xmlName); 
    }
} 
