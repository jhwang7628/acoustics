#include <iostream>
#include <string>
#include <parser/RigidsimParseConfig.h>
#include <rigid/RigidsimConfigData.h>

int main()
{
    std::cout << "test_Rigid\n";
    const std::string configFile("/home/jui-hsien/code/acoustics/work/plate_drop_test/default.cfg"); 
    RigidsimParseConfig parser(configFile);
    std::cout << "test_Rigid 1\n";
    parser.Parse(); 
    std::cout << "test_Rigid 2\n";
    auto data = parser.GetParsedData(); 
    std::cout << "test_Rigid 3\n";
    std::cout << data->simulation_step_size << " " << data->simulation_time_len << std::endl;
    return 0; 
}
