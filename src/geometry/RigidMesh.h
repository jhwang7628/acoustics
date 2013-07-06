//////////////////////////////////////////////////////////////////////
// RigidMesh.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef RIGID_MESH_H
#define RIGID_MESH_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/Matrix3.hpp>
#include <linearalgebra/VECTOR.h>
#include <linearalgebra/Vector3.hpp>

#include <geometry/Point3.hpp>
#include <geometry/TriangleMesh.hpp>

#include <TYPES.h>

//////////////////////////////////////////////////////////////////////
// RigidMesh class
//
// Stores a triangle mesh, and relevant rigid body data
//////////////////////////////////////////////////////////////////////
class RigidMesh {
    public:
        RigidMesh( const TriangleMesh<REAL> &mesh, const std::string &rigidFilePrefix,
                   REAL density );

        // Destructor
        virtual ~RigidMesh();

        // Figures out an inverse effective mass for the given
        // position on the surface of the mesh, and a force applied in the
        // given direction
        REAL inverseEffectiveMass( const Point3d &x,
                                   const Vector3d &forceDirection ) const;

        // Returns acceleration directional magnitudes, given a force
        // applied in the given direction to the given point
        void accelerationCoefficients( const Point3d &x, const Vector3d &force,
                                       Vector3d &translationalAcceleration,
                                       Vector3d &rotationalAcceleration ) const;

        REAL mass() const
        {
            return _mass;
        }

        REAL youngsModulus() const
        {
            return _youngsModulus;
        }

        REAL poissonRatio() const
        {
            return _poissonRatio;
        }

        const Point3d &centerOfMass() const
        {
            return _centerOfMass;
        }

        const Matrix3d &inertiaAxes() const
        {
            return _inertiaAxes;
        }

        REAL ellipseScaleX() const
        {
            return _xScale;
        }

        REAL ellipseScaleY() const
        {
            return _yScale;
        }

        REAL ellipseScale() const
        {
            return _matchingScale;
        }

        const TriangleMesh<REAL> &mesh() const
        {
            return _mesh;
        }

        // Fits this mesh to an ellipsoid in which the z-axis has unit length
        void fitEllipsoid( ObjectMeasure measure = VOLUME );

        // Resets the fit ellipsoid for this mesh and determines a new scaling
        void setEllipsoid( REAL xScale, REAL yScale,
                           ObjectMeasure measure = VOLUME );

    protected:

    private:
        const TriangleMesh<REAL>    &_mesh;

        // Rigid body quantities
        REAL                         _mass;
        REAL                         _density;
        Point3d                      _centerOfMass;
        Matrix3d                     _inertia;
        Matrix3d                     _inertiaInverse;

        Matrix3d                     _inertiaAxes;
        Vector3d                     _inertiaMoments;

        // Material properties
        REAL                         _youngsModulus;
        REAL                         _poissonRatio;

        // Ellipsoid properties
        REAL                         _xScale;
        REAL                         _yScale;

        REAL                         _matchingScale;

};

#endif
