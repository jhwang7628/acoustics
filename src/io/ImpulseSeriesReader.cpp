#include <fstream>
#include <map>
#include <io/ImpulseSeriesReader.h>
#include <parser/RigidsimParseConfig.h>

//##############################################################################
//##############################################################################
void ImpulseSeriesReader::
LoadRigidsimConfig(ImpulseSeriesObjectPtr object)
{
    assert(IO::ExistFile(_rigidsimConfigFile));
    RigidsimParseConfig parser(_rigidsimConfigFile);
    parser.Parse();
    object->SetRigidsimConfigData(parser.GetParsedData());
}

//##############################################################################
// Load the impulses from file for object identified. The format is the following:
//
//  <timestamp>  <object id(0-based)>  <vertex ID>  <relative speed> <impulse.x>
//  <impulse.y> <impulse.z>  <T/S> <C/P>
//
//  There are two types of impulse ("C/P"): contraint impulse (C) that records
//  impulses when object hit the ground (a constraint), and pair impulse (P)
//  that records impulses between two interacting objects. Impulses are defined
//  in the object frame of reference.
//
//  The letter "T/S" is used to designate the different object in the pair
//  impulse case. "S" is used for the rigid body of class TRigidBody ba in
//  RigidObjImpRecorder::record_inter_obj_impulse, and "T" is used for the
//  rigid body of class TRigidBody bb in the same function.
//
//  Note:
//   The file is written by the class
//   io/RigidObjImpRecorder::record_constraint_impulse to object m_modalImpulseOut
//
//##############################################################################
void ImpulseSeriesReader::
LoadImpulses(const int &loadObjectID, ImpulseSeriesObjectPtr object, std::shared_ptr<FDTD_Objects> objects)
{
    std::ifstream inFile(_impulseFile);
    if (!inFile)
    {
        std::cout << "**WARNING** Cannot open impulse file: " + _impulseFile
                  << std::endl << "Skipping..\n";
        return;
    }

    RigidSoundObjectPtr soundObject =
        std::dynamic_pointer_cast<FDTD_RigidSoundObject>(objects->GetPtr(loadObjectID));
    int objectID=-1, objectID_old=-1, count=0;
    char pairOrder, impulseType;
    ImpulseSeriesObject::ImpactRecord buffer, buffer_old;
    while(inFile >> buffer.timestamp >> objectID >> buffer.appliedVertex >> buffer.contactSpeed
                 >> buffer.impactVector.x >> buffer.impactVector.y >> buffer.impactVector.z
                 >> pairOrder >> impulseType)
    {
        if (objectID == loadObjectID)
        {
            if (impulseType == 'C')
            {
                // estimate contact time scale between this object and
                // constraint, and added to impulse train
                buffer.supportLength = soundObject->EstimateContactTimeScale(buffer.appliedVertex, buffer.contactSpeed, buffer.impactVector);
                buffer.impactPosition = Vector3d(soundObject->GetMeshPtr()->vertex(buffer.appliedVertex));
                buffer.gamma = M_PI*buffer.impactVector.length()/2.0/buffer.supportLength;
                object->AddImpulse(buffer);  // assumes point impulse
            }
            else if (impulseType == 'P')
            {
                if (pairOrder == 'T') // old buffer contains the pair 'S'
                {
                    assert(count > 0); // make sure buffer exists
                    RigidSoundObjectPtr pairObject =
                        std::dynamic_pointer_cast<FDTD_RigidSoundObject>(objects->GetPtr(objectID_old));
                    buffer.supportLength = soundObject->EstimateContactTimeScale(pairObject, buffer.appliedVertex, buffer_old.appliedVertex, buffer.contactSpeed, buffer.impactVector);
                    buffer.impactPosition = Vector3d(soundObject->GetMeshPtr()->vertex(buffer.appliedVertex));
                    buffer.gamma = M_PI*buffer.impactVector.length()/2.0/buffer.supportLength;
                    object->AddImpulse(buffer); // add the one that's in buffer
                }
                else if (pairOrder == 'S') // need to read one more line to get pair
                {
                    buffer_old = buffer;
                    objectID_old = objectID;
                    inFile >> buffer.timestamp >> objectID >> buffer.appliedVertex >> buffer.contactSpeed
                           >> buffer.impactVector.x >> buffer.impactVector.y >> buffer.impactVector.z
                           >> pairOrder >> impulseType;
                    RigidSoundObjectPtr pairObject =
                        std::dynamic_pointer_cast<FDTD_RigidSoundObject>(objects->GetPtr(objectID));
                    buffer_old.supportLength = soundObject->EstimateContactTimeScale(pairObject, buffer_old.appliedVertex, buffer.appliedVertex, buffer_old.contactSpeed, buffer_old.impactVector);
                    buffer_old.impactPosition = Vector3d(soundObject->GetMeshPtr()->vertex(buffer_old.appliedVertex));
                    buffer_old.gamma = M_PI*buffer_old.impactVector.length()/2.0/buffer_old.supportLength;
                    object->AddImpulse(buffer_old); // add the one that's in buffer_old (one labeled 'S')
                }
            }
        }
        buffer_old = buffer;
        objectID_old = objectID;
        count ++;
    }
    object->Filter();
    LoadRigidsimConfig(object);
}

//##############################################################################
// This function loads impulses but cannot compute contact time scales and thus
// all impulses are assumed to have zero duration (point impulse). Suitable if
// only modal vibration is needed.
//##############################################################################
void ImpulseSeriesReader::
LoadImpulses(const int &loadObjectID, ImpulseSeriesObjectPtr object)
{
    std::ifstream inFile(_impulseFile);
    if (!inFile)
        throw std::runtime_error("**ERROR** Cannot open impulse file: " + _impulseFile);

    int objectID;
    char pairOrder, impulseType;
    ImpulseSeriesObject::ImpactRecord buffer;
    while(inFile >> buffer.timestamp >> objectID >> buffer.appliedVertex >> buffer.contactSpeed
                 >> buffer.impactVector.x >> buffer.impactVector.y >> buffer.impactVector.z
                 >> pairOrder >> impulseType)
    {
        if (objectID == loadObjectID)
        {
            if (impulseType == 'C' || impulseType == 'P')
            {
                // note: here I assume impulse has 0 support (default)
                object->AddImpulse(buffer);
            }
        }
    }
    object->Filter();
    LoadRigidsimConfig(object);
}

