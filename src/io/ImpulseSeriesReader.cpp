#include <fstream>
#include <map>
#include <io/ImpulseSeriesReader.h> 

//##############################################################################
// Load the impulses from file. The format is the following: 
//
//  <timestamp>  <object id(0-based)>  <vertex ID>  <impulse.x>  <impulse.y>
//  <impulse.z>  <T/S>
//  
//  The last letter "T/S" indicates what kind of vertex ID is used. "T" means
//  the vertex ID is the id in tetrahedron mesh, whereas "S" means the vertex 
//  ID is the id in surface triangle mesh.
//
//  (comments copied from sndgen/RigidSoundObj.hpp)
//
//  Note: 
//   The file is written by the class
//   io/RigidObjImpRecorder::record_constraint_impulse to object m_modalImpulseOut
//
//##############################################################################
void ImpulseSeriesReader::
LoadImpulses(std::vector<ImpulseSeriesObjectPtr> &objects)
{
    // first check if file exists 
    std::ifstream inFile(_impulseFile); 
    if (!inFile) 
        throw std::runtime_error("**ERROR** Cannot open impulse file: " + _impulseFile);

    // now go through the first couple of lines of the file to see if its
    // format is correct, and check the dimension of the passed-in objects
    int count = 0; 
    std::map<int, int> objectIDMap;  // key: object id from file. value: relative position within first timestamp.
    
    double timestamp, firstTimestamp;
    int objectID, vertexID; 
    Vector3d impulse; 
    char impulseType, endLine; 
    while(inFile >> timestamp >> objectID  >> vertexID
                 >> impulse.x >> impulse.y >> impulse.z
                 >> impulseType >> endLine && (firstTimestamp == timestamp || count==0)) 
    {
        if (impulseType != 'S') 
            throw std::runtime_error("**ERROR** Impulse vertex is prescribed on tetrahedron mesh, this case is not handled."); 
        if (endLine != 'C')
            throw std::runtime_error("**ERROR** End of line character is not 'C', check file format."); 
        if (count == 0) 
            firstTimestamp = timestamp; 

        // check if the id exists in the map
        auto search = objectIDMap.find(objectID); 
        if (search != objectIDMap.end())
            throw std::runtime_error("**ERROR** Within the first timestamp, there are two impulses registered with the same id. This is not allowed."); 
        else 
            objectIDMap[objectID] = count; 
        count ++; 
    }
    inFile.close(); 

    // pass the test, actually read the file and add impulses to objects.
    const bool reinitializeObjects = (objects.size()!=(size_t)count ? true : false); 
    if (reinitializeObjects)
    {
        std::cout << "Initializing objects that stores the impulses.\n";
        objects.resize(count); 
    }
    inFile.open(_impulseFile);
    while(inFile >> timestamp >> objectID  >> vertexID
                 >> impulse.x >> impulse.y >> impulse.z
                 >> impulseType >> endLine) 
    {
        ImpulseSeriesObjectPtr &object = objects.at(objectID); 
        if (reinitializeObjects) 
            object = std::make_shared<ImpulseSeriesObject>(objectID, nullptr);
        object->AddImpulse(timestamp, vertexID, impulse); 
    }
}
