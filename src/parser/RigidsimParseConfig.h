#ifndef RIGIDSIM_PARSE_CONFIG_H 
#define RIGIDSIM_PARSE_CONFIG_H 
#include <memory>
#include <parser/ParseConfig.h> 
#include <rigid/RigidsimConfigData.h> 

//##############################################################################
// This class helps parsing the config file for tools/rigidsim
//##############################################################################
class RigidsimParseConfig : public ParseConfig
{
    private: 
        std::shared_ptr<RigidsimConfigData> _parsedData; 

    public: 
        RigidsimParseConfig() 
            : ParseConfig()
        {}
        RigidsimParseConfig(const std::string &configFile) 
            : ParseConfig(configFile)
        {}

        inline std::shared_ptr<RigidsimConfigData> GetParsedData(){return _parsedData;}
        virtual void Parse(); 

};

#endif
