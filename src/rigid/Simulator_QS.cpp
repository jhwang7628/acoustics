#include "Simulator_QS.h"
#include "utils/print_msg.h"

using namespace std;

/*
 * this method should be called each time before the advance method is called
 */
void RigidSimulator_QS::begin_advance()
{
    //// save the velocity in the last timestep
    const map<TRigidBody*, RigidStateRec>::iterator end = stateRec_.end();
    for(map<TRigidBody*, RigidStateRec>::iterator it = stateRec_.begin();
            it != end;++ it)
    {
        it->second.lastVel = it->first->velocity();
        it->second.lastAngVel = it->first->angular_velocity();

        it->second.lastPos = it->first->mass_center();
        it->second.lastRot = it->first->rotation();
    }
}

/*
 * this method should be called each time after the advance method is called
 *
 * it restores the old position and velocity of the quasi-static objects
 */
void RigidSimulator_QS::end_advance()
{
#ifdef WITH_QUASI_STATIC_OBJS
    const map<TRigidBody*, RigidStateRec>::iterator end = stateRec_.end();
    for(map<TRigidBody*, RigidStateRec>::iterator it = stateRec_.begin();
            it != end;++ it)
        if ( it->second.isStatic )
        {
            it->first->set_position(it->second.lastPos, it->second.lastRot);
            it->first->set_velocity(it->second.lastVel, it->second.lastAngVel);
        }
#endif
    ////// damp the state
    //if ( time() > 0.003 )
    //{
    //    cerr << "!!!!!! Dump state to state.dump" << endl;
    //    ofstream fout("state.dump");
    //    fout << setprecision(20);
    //    for(size_t i = 0;i < m_rigidObjs.size();++ i)
    //    {
    //        const Point3d mc = m_rigidObjs[i]->mass_center();
    //        const Quat4d rot = m_rigidObjs[i]->rotation();

    //        fout << m_rigidObjs[i]->id() << ' '
    //             << mc.x << ' ' << mc.y << ' ' << mc.z << ' '
    //             << rot.w << ' ' << rot.v.x << ' ' << rot.v.y << ' ' << rot.v.z << endl;
    //    }
    //    fout.close();
    //}
}

void RigidSimulator_QS::set_rigid_state(
        int id, const Tuple3d& pos, const Quat4d& rot)
{
    if ( id >= m_rigidObjs.size() ) return;
    if ( m_rigidObjs[id]->id() != id )
    {
        PRINT_ERROR("Inconsistent ID from dump file\n");
        exit(1);
    }
    m_rigidObjs[id]->set_position(pos, rot);
}

void RigidSimulator_QS::set_debris_state(TRigidBody* bd, 
        TRigidBody* parent, bool isstatic) 
{
    assert(!stateRec_.count(parent));

#ifdef WITH_QUASI_STATIC_OBJS
    // set the current position
    if ( isstatic )
    {
        const RigidStateRec& rec = stateRec_[parent];
        Point3<REAL> pos = rec.lastPos + rec.lastRot.rotate(
                bd->initial_mass_center() - parent->initial_mass_center());
        bd->set_position(pos, rec.lastRot);
    }
    else
#endif
    {
        Point3<REAL> pos = parent->current_position(
                bd->initial_mass_center());
        bd->set_position(pos, parent->rotation());
    }
}

void RigidSimulator_QS::set_last_velocity(TRigidBody* bd, TRigidBody* parent)
{
    const RigidStateRec& rec = stateRec_[parent];
    Vector3<REAL> lastvec = rec.lastRot.rotate(
                bd->initial_mass_center() - parent->initial_mass_center());
    bd->set_velocity(rec.lastVel + rec.lastAngVel.crossProduct(lastvec),
                     rec.lastAngVel);
}

