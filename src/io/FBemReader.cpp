#include <io/FBemReader.h> 

//##############################################################################
//##############################################################################
bool FBemReader::
CheckFBemInputAgainstMesh(std::shared_ptr<TriangleMesh<REAL> > &mesh, const std::string &solutionFile)
{
    std::cout << "check input\n";
}

//##############################################################################
//##############################################################################
bool FBemReader::
ReadFBemOutputToInfo(std::shared_ptr<BEMSolutionMode> &solution, const std::string &solutionFile)
{
    std::cout << "check output\n";
}
