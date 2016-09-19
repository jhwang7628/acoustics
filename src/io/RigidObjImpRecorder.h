#ifndef IO_RIGID_IMPULSE_RECORDER_H
#   define IO_RIGID_IMPULSE_RECORDER_H

#include <fstream>
#include "geometry/FixVtxTetMesh.hpp"
#include "rigid/LSCollisionRigidBody.hpp"
#include "rigid/CollisionRec.hpp"

/*
 * Record impulses applied on objects into file. These impulses 
 * will be used to drive the vibration later when generating the 
 * sound.
 *
 * The recording in the file is gonna to be
 * <time>    <obj id>    <applied vtx id>    <impulse>
 * ...
 *
 * NOTE that the "applied vtx id" is the ID on surface triangle mesh
 *
 */
class RigidObjImpRecorder
{
    public:
        typedef FixVtxTetMesh<REAL>                     TMesh;
        typedef LSCollisionRigidBody<REAL, TMesh>       TRigidBody;

    public:
        ~RigidObjImpRecorder()
        {
            if ( m_fout.is_open() ) {
                m_fout.close();
            }

            if ( m_modalImpulseOut.is_open() ) {
                m_modalImpulseOut.close();
            }
        }

        void init(const char* file, const char *modalFile, int precision = 16);
        /*
         * record the impulse, applied on the vertex (id = cRec.vtxId)
         * of the rigid object (body).
         */
        void record_impulse(REAL ts, const Vector3<REAL>& imp, 
                int vtxId, const TRigidBody* body, bool surfVtx);
        void record_impulse(REAL ts, const Vector3<REAL>& imp, 
                const Point3<REAL>& pt, const TRigidBody* body);

        void record_constraint_impulse(REAL ts, 
                const Vector3<REAL>& imp, const CollisionRec<REAL>& cRec,
                const TRigidBody* body);

        void record_inter_obj_impulse(REAL ts, 
                const Vector3<REAL>& imp, const CollisionRec<REAL>& cRec,
                const TRigidBody* ba, const TRigidBody* bb);

    private:
        std::ofstream       m_fout;

        std::ofstream       m_modalImpulseOut;
};

#endif
