#include <modal_model/ModalAnalysis.h> 

//##############################################################################
//##############################################################################
void ModalAnalysis::
BuildModalModelsFromFile()
{
    std::cout << "Load configuration file: " << _configFile << std::endl;
    libconfig::Config config; 
    try
    {
        config.readFile(_configFile.c_str()); 
        libconfig::Setting &sndobjSetting = config.lookup("sndobj"); 
        const int N_sndobj = sndobjSetting.getLength(); 
        for (int obj=0; obj<N_sndobj; ++obj)
        {
            RigidModalPtr rigidModal = std::make_shared<RigidModal>(sndobjSetting[obj]);
            _rigidModals.push_back(rigidModal); 
        }
        std::cout << "Number of modal models built: " << _rigidModals.size() << std::endl;
    }
    catch (const libconfig::SettingException& e)
    {
        std::cerr << "**ERROR** Error reading configuration file at " << e.getPath() << ": " << e.what() << "\n";
        return; 
    }
    catch (libconfig::ParseException& e)
    {
        std::cerr << "**ERROR** Error parsing configuration file : " << e.getError() << " " << e.what() << " at Line " << e.getLine() << "\n";
        return; 
    }
}
