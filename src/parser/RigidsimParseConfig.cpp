#include <string>
#include <parser/RigidsimParseConfig.h> 

//##############################################################################
//##############################################################################
void RigidsimParseConfig::
Parse()
{
    libconfig::Config config;
    try
    {
        config.readFile(_configFile.c_str()); 
        _parsedData = std::make_shared<RigidsimConfigData>();
        LOOKUP_VALUE_GUARD(config, "simulation.step_size", _parsedData->simulation_step_size); 
        LOOKUP_VALUE_GUARD(config, "simulation.time_len", _parsedData->simulation_time_len); 
    }
    catch (const libconfig::SettingException &e)
    {
        fprintf(stderr, "Error occured when reading configure file at %s: %s\n",
                e.getPath(), e.what());
    }
    catch (const libconfig::ParseException &e)
    {
        fprintf(stderr, 
                "ParseException: %s %s at Line %d\n",
                e.getError(), e.what(), e.getLine());
    }

    if (!_parsedData)
        throw std::runtime_error("**ERROR** Cannot parse rigidsim config: " + _configFile); 
}
