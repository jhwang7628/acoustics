#ifndef DIFF_DEFINE
/******************************************************************************
 *  File: RigidObjDispRecorder.h
 *
 *  This file is part of isostuffer
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#endif /* ! DIFF_DEFINE */
#ifndef IO_RIGID_DISP_RECORDER_H
#define IO_RIGID_DISP_RECORDER_H

#include <fstream>
#include "config.h"
#include "geometry/FixVtxTetMesh.hpp"
#include "rigid/LSCollisionRigidBody.hpp"

/*
 * Record the displacement (translation and rotation into files)
 *
 * <d_ts> 
 * <i_id1>   # first rigid object 
 * <d_translate.x> <d_translate.y> <d_translate.z> 
 * <d_rot.w> <d_rot.x> <d_rot.y> <d_rot.z>
 * <i_id2>   # 2nd rigid obj
 * <d_translate.x> <d_translate.y> <d_translate.z> 
 * <d_rot.w> <d_rot.x> <d_rot.y> <d_rot.z>
 * <-1>      # label the ending of this timestep
 * ...
 * <d_ts>
 * <i_id1>   # first rigid object ...
 *
 *
 * Record the velocity (translation and rotation into files)
 *
 * <d_ts> 
 * <i_id1>   # first rigid object 
 * <d_velocity.x> <d_velocity.y> <d_velocity.z> 
 * <d_angular_velocity.x> <d_angular_velocity.y> <d_angular_velocity.z> 
 * <i_id2>   # 2nd rigid obj
 * <d_velocity.x> <d_velocity.y> <d_velocity.z> 
 * <d_angular_velocity.x> <d_angular_velocity.y> <d_angular_velocity.z> 
 * <-1>      # label the ending of this timestep
 * ...
 * <d_ts>
 * <i_id1>   # first rigid object ...
 *
 *
 * Record the velocity (translation and rotation into files)
 *
 * <d_ts> 
 * <i_id1>   # first rigid object 
 * <d_velocity.x> <d_velocity.y> <d_velocity.z> 
 * <d_angular_velocity.x> <d_angular_velocity.y> <d_angular_velocity.z> 
 * <i_id2>   # 2nd rigid obj
 * <d_velocity.x> <d_velocity.y> <d_velocity.z> 
 * <d_angular_velocity.x> <d_angular_velocity.y> <d_angular_velocity.z> 
 * <-1>      # label the ending of this timestep
 * ...
 * <d_ts>
 * <i_id1>   # first rigid object ...
 */
class RigidObjDispRecorder
{
    public:
        typedef FixVtxTetMesh<REAL>                     TMesh;
        typedef LSCollisionRigidBody<REAL, TMesh>       TRigidBody;
        enum RecordFlag{ DISP_ONLY=0, DISP_VELO_ACCE=1 }; 

        void init(const char* file)
        {
            if ( m_fout.is_open() ) {
                m_fout.close();
            }
            m_fout.open(file, std::ios::binary);
            m_record_flag = DISP_ONLY; 
        }

        void init(const char* file_disp, const char* file_velo, const char* file_acce)
        {
            if ( m_fout.is_open() ) {
                m_fout.close();
            }
            m_fout.open(file_disp, std::ios::binary);

            if ( m_fout_velo.is_open() ) {
                m_fout_velo.close();
            }
            m_fout_velo.open(file_velo, std::ios::binary);

            if ( m_fout_acce.is_open() ) {
                m_fout_acce.close();
            }
            m_fout_acce.open(file_acce, std::ios::binary);
            m_record_flag = DISP_VELO_ACCE; 
        }

        /*
         * Begin recording the displacement (both translation and rotation)
         * for the given time moment
        */
        void begin_record(double ts)
        {
            m_fout.write((const char*)&ts, sizeof(double));
            if (m_record_flag == DISP_VELO_ACCE) 
            {
                m_fout_velo.write((const char*)&ts, sizeof(double));
                m_fout_acce.write((const char*)&ts, sizeof(double));
            }
        }

        /*
         * write the displacement of the given rigid body into file
         */
        void record_displacement(const TRigidBody* body)
        {
            const int id = body->id();
            m_fout.write((const char*)&id, sizeof(int));

            const Point3<REAL>& c = body->mass_center();
            m_fout.write((const char*)&c, sizeof(Point3<REAL>));
            
            const Quaternion<REAL>& r = body->rotation();
            m_fout.write((const char*)&r, sizeof(Quaternion<REAL>));
        }

        /*
         * write the displacement, velocity, and acceleration of the given rigid body into file
         */
        void record_rigid_kinematics(const TRigidBody* body) 
        {
            assert(m_fout && m_fout_velo && m_fout_acce); // files have to be defined 

            // write displacement
            record_displacement(body); 

            // write velocity
            const int id = body->id();
            m_fout_velo.write((const char*)&id, sizeof(int));

            const Vector3<REAL>& v = body->velocity();
            m_fout_velo.write((const char*)&v, sizeof(Vector3<REAL>));
            
            const Vector3<REAL>& o = body->angular_velocity();
            m_fout_velo.write((const char*)&o, sizeof(Vector3<REAL>));

            // write accleration
            m_fout_acce.write((const char*)&id, sizeof(int));

            const Vector3<REAL>& a = body->acceleration();
            m_fout_acce.write((const char*)&a, sizeof(Vector3<REAL>));
            
            const Vector3<REAL>& k = body->angular_acceleration();
            m_fout_acce.write((const char*)&k, sizeof(Vector3<REAL>));
        }

        void end_record()
        {
            const int END = -1;
            m_fout.write((const char*)&END, sizeof(int));
            // flush data after each time step
            m_fout.flush();
            if (m_record_flag == DISP_VELO_ACCE) 
            {
                m_fout_velo.write((const char*)&END, sizeof(int));
                m_fout_velo.flush();
                m_fout_acce.write((const char*)&END, sizeof(int));
                m_fout_acce.flush();
            }
        }

        ~RigidObjDispRecorder()
        {
            if ( m_fout.is_open() )
                m_fout.close();
            if (m_record_flag == DISP_VELO_ACCE && m_fout_velo.is_open())
                m_fout_velo.close();
            if (m_record_flag == DISP_VELO_ACCE && m_fout_acce.is_open())
                m_fout_velo.close();
        }

    private:
        RecordFlag          m_record_flag; 
        std::ofstream       m_fout;      // displacement file
        std::ofstream       m_fout_velo; // velocity file
        std::ofstream       m_fout_acce; // acceleration file
};

#endif
