#ifndef PARSE_CONFIG_H
#define PARSE_CONFIG_H
#include <libconfig.h++>
#include <stdexcept>
//##############################################################################
// Macros (only for this and inherited class
//##############################################################################
#ifndef LOOKUP_VALUE_GUARD
#define LOOKUP_VALUE_GUARD(config, varName, varReference) \
    if (!config.lookupValue(varName, varReference)) \
            throw std::runtime_error("**ERROR** Cannot find "+std::string(varName)+" in config file"); 
#endif


//##############################################################################
// Helper class for libConfig
//##############################################################################
class ParseConfig
{
    protected:
        std::string _configFile; 

    public: 
        ParseConfig()
        {}
        ParseConfig(const std::string &configFile)
            : _configFile(configFile)
        {}

        inline std::string GetConfigFile(){return _configFile;}
        inline const std::string GetConfigFile() const {return _configFile;}
        inline void SetConfigFile(const std::string &configFile){_configFile = configFile;} 
        virtual void Parse() = 0; 
};

#endif
