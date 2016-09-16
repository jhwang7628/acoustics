#ifndef RIGID_SIMULATOR_HPP
#   define RIGID_SIMULATOR_HPP

#include "constants.h"
#include "geometry/FixVtxTetMesh.hpp"
#include "LSCollisionRigidBody.hpp"
#include "CollisionRec.hpp"
#include "ContactGraph.h"

#ifdef USE_RECORDER
#include "io/RigidObjImpRecorder.h"
#include "io/RigidObjDispRecorder.h"
#endif

class RigidBodySimulator
{
    public:
        typedef FixVtxTetMesh<REAL>                     TMesh;
        typedef LSCollisionRigidBody<REAL, TMesh>       TRigidBody;
        typedef TRigidBody::TCollProc                   TCollProc;
        typedef CollisionConstraint<REAL, TCollProc>    TConstraint;
#ifdef USE_RECORDER
        typedef RigidObjImpRecorder                     TImpRecorder;
#endif
        // this struct stores if the rigid body motion is restricted on a plane
        struct MotionProjection
        {
            bool            useProjection; 
            Vector3<REAL>   projectionNormal; // normal of the plane the motion is allowed on
            MotionProjection() : useProjection(false){}
        };

    public:
        RigidBodySimulator();
#ifdef USE_RECORDER
        RigidBodySimulator(const char* recFile, const char *modalRecFile,
                           const char* dispFile);
#endif
        ~RigidBodySimulator();

        /*
         * NOTE: this method should be called each time when some old objects
         *       are removed or some new ones are added.
         */
        void init();

        void advance(REAL dt);

        const std::vector<TRigidBody*>& rigid_bodies() const
        {   return m_rigidObjs; }

        inline void add_rigid_body(TRigidBody* bd)
        {
          if ( bd->mesh()->num_fixed_vertices() > 0 )
          {
            bd->set_fixed( true );
          }
          m_rigidObjs.push_back(bd);
        }
        inline void add_constraint(TConstraint* ct)
        {  m_constraints.push_back(ct); }

        TRigidBody* rigid_body(int i) 
        {  return m_rigidObjs[i]; }

        /*
         * NOTE: after calling this method, make sure init() is called before the 
         *       next time step.
         */
        inline void remove_rigid_body(TRigidBody* bd)
        {
            std::vector<TRigidBody*>::iterator end = m_rigidObjs.end();
            m_rigidObjs.erase(std::remove(m_rigidObjs.begin(), end, bd), end);
        }

        inline void assign_motion_projection(const Vector3<REAL> &projectionNormal)
        {
            m_motionProjection.projectionNormal = projectionNormal; 
            m_motionProjection.useProjection = true;
        }

#ifdef USE_RECORDER
        TImpRecorder& impulse_recorder()
        {  return m_impRec; }

        TImpRecorder* impulse_recorder_ptr()
        {  return &m_impRec; }
#endif

        REAL time() const
        {   return m_ts; }

        void set_time(REAL ts)
        {   m_ts = ts;  }

    private:
        void collision_matrix(const TRigidBody* b, 
                const Vector3<REAL>& r, Matrix3<REAL>& K);

        /* ========== For Collisions ========== */
        void resolve_collision(REAL dt);
        int  collision_pair_resp(TRigidBody* b1, TRigidBody* b2, REAL dt);
        /* apply collision constraints on rigid body */
        int  collision_constraint_resp(TRigidBody* b, REAL dt);
        /* apply impulse to the rigid bodies */
        void apply_collision_impulse(TRigidBody* b1, TRigidBody* b2, 
                    CollisionRec<REAL>& cRec);
        void apply_collision_impulse(TRigidBody* b, 
                    CollisionRec<REAL>& cRec);

        /* ========== For Contacts ========== */
        void resolve_contact(REAL dt);
        int  contact_pair_resp(TRigidBody* b1, TRigidBody* b2, REAL dt, int numIt);
        /* apply collision constraints on rigid body */
        int  contact_constraint_resp(TRigidBody* b, REAL dt, int numIt);

        void apply_contact_impulse(TRigidBody* b1, TRigidBody* b2, 
                    const CollisionRec<REAL>& cRec, const REAL restC);
        void apply_contact_impulse(TRigidBody* b, 
                    const CollisionRec<REAL>& cRec, const REAL restC);

    protected:
        std::vector<TRigidBody*>    m_rigidObjs;

    private:
        std::vector<int>            m_shuffledObjIds;
        std::vector<TConstraint*>   m_constraints;  // the constraints like walls
        ContactGraph                m_contactGraph;
        REAL                        m_ts;
        MotionProjection            m_motionProjection;

#ifdef USE_RECORDER
        TImpRecorder                m_impRec;       // impulse recorder
        RigidObjDispRecorder        m_dispRec;      // displacement recorder
#endif
};

#endif

