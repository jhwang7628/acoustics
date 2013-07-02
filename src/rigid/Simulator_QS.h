#ifndef RIGID_SIMULATOR_QUASI_STATIC_OBJECTS_H
#   define RIGID_SIMULATOR_QUASI_STATIC_OBJECTS_H

#include "Simulator.h"
#include "utils/print_msg.h"

/*!
 * Rigid body simulator with Quasi-static objects
 *
 * When no WITH_QUASI_STATIC_OBJS is defined, this simulator just remember the velocity
 * of each object at last timestep, this information is used to compute the velocity
 * of fractured objects by fracture simulator
 *
 * When WITH_QUASI_STATIC_OBJS is defined, this class also provides the methods to 
 * support the quasi-static objects. Quasi-static object is the object which is fixed
 * (static) but when simulate them, we treat them as free object, but fix their 
 * positions at the end of each time step
 */
class RigidSimulator_QS : public RigidBodySimulator
{
    public:
        typedef RigidBodySimulator::TRigidBody      TRigidBody;

        struct RigidStateRec
        {
            Vector3<REAL>       lastVel;        // velocity at last time step
            Vector3<REAL>       lastAngVel;     // angular velocity at last time step

            Point3<REAL>        lastPos;
            Quaternion<REAL>    lastRot;
#ifdef WITH_QUASI_STATIC_OBJS
            bool                isStatic;
#endif
        };

    public:
        RigidSimulator_QS()
        { }
#ifdef USE_RECORDER
        RigidSimulator_QS(const char* recFile, const char* modalRecFile,
                          const char* dispFile):
                RigidBodySimulator(recFile, modalRecFile, dispFile)
        { }
#endif

        /*!
         * Note that, by default, the objects with fixed vertices are always
         * quasi-static objects
         */
        void add_rigid_body(TRigidBody* bd)
        {
            assert(!stateRec_.count(bd));

#ifdef WITH_QUASI_STATIC_OBJS
            stateRec_[bd].isStatic = bd->mesh()->num_fixed_vertices() > 0;
            if ( stateRec_[bd].isStatic ) PRINT_MSG("static object is added\n");
#endif
            stateRec_[bd].lastVel.set(0, 0, 0);
            stateRec_[bd].lastAngVel.set(0, 0, 0);
            RigidBodySimulator::add_rigid_body(bd);
        }

#ifdef WITH_QUASI_STATIC_OBJS
        /*
         * Add an object with an indication about that if the object is 
         * quasi-static or not
         */
        void add_rigid_body(TRigidBody* bd, bool isquasi)
        {
            assert(!stateRec_.count(bd));

            stateRec_[bd].isStatic = isquasi;
            stateRec_[bd].lastVel.set(0, 0, 0);
            stateRec_[bd].lastAngVel.set(0, 0, 0);
            RigidBodySimulator::add_rigid_body(bd);
        }

        /*
         * given a new point in the object, bd, inital configuration space, return 
         * the point in the last configuration space of that object
         */
        inline Point3<REAL> last_position(TRigidBody* bd, const Point3<REAL>& pt) 
        {
            assert(!stateRec_.count(bd));

            const RigidStateRec& rec = stateRec_[bd];
            return rec.lastPos + rec.lastRot.rotate(pt - bd->initial_mass_center());
        }
#endif

        /*
         * Add a new rigid body which is debris of the rigid object, parent.
         * So this method adds the new rigid body and set its current position
         * and rotation as its parent's last position and rotation
         */
        void set_debris_state(TRigidBody* bd, TRigidBody* parent, bool isstatic);
        void set_last_velocity(TRigidBody* bd, TRigidBody* parent);

        std::pair< const Vector3<REAL>*, const Vector3<REAL>* >
        last_velocity(TRigidBody* pbody)
        {
            assert(pbody && stateRec_.count(pbody));
            const RigidStateRec& rec = stateRec_[pbody];
            return std::make_pair(&(rec.lastVel), &(rec.lastAngVel));
        }

        void begin_advance();
        void end_advance();

        void set_rigid_state(int id, const Tuple3d& pos, const Quat4d& rot);

        /*
         * NOTE: after calling this method, make sure init() is called before the 
         *       next time step.
         */
        void remove_rigid_body(TRigidBody* bd)
        {
            stateRec_.erase(bd);
            RigidBodySimulator::remove_rigid_body(bd);
        }

    private:
        std::map<TRigidBody*, RigidStateRec>    stateRec_;
};

#endif
