#include <string>
#include <iostream>
#include <utils/STL_Wrapper.h> 
#include <deformable/ModeData.h> 

void TestRigidBodySim()
{
    const std::string testDisplacement("/home/jui-hsien/code/acoustics/work/plate_drop_test/displace.bin");
}

void TestIO()
{
    const std::string filename("/home/jui-hsien/code/acoustics/work/ball_0p1/small_ball_r_0p10.modes"); 
    ModeData modeData; 
    modeData.read(filename.c_str()); 
    std::cout << modeData << std::endl;
    for (int mm=0; mm<modeData.numModes(); ++mm) 
    {
        std::cout << "Mode " << mm << ": " << modeData._omegaSquared.at(mm) << std::endl;
        std::cout << " shape: "; 
        STL_Wrapper::PrintVectorContent(std::cout, modeData.mode(mm), 8); 
    }
}

int main() 
{
    std::cout << "Unit Test: Modal Sound\n"; 
    TestIO(); 
    return 0; 
}
