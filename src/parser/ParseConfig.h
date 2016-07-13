#ifndef PARSE_CONFIG
#define PARSE_CONFIG
#include <libconfig.h++>

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
