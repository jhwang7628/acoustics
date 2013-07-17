//////////////////////////////////////////////////////////////////////
// GTS_TriMesh.h: Interface for the GTS_TriMesh class
//
//////////////////////////////////////////////////////////////////////

#ifndef GTS_TRI_MESH_H
#define GTS_TRI_MESH_H

#include <TYPES.h>

#include "TriangleMesh.hpp"

#include <linearalgebra/Vector3.hpp>

#include <gts.h>

//////////////////////////////////////////////////////////////////////
// GTS_TriMesh class
//
// GTS representation for our tri mesh class
//////////////////////////////////////////////////////////////////////
class GTS_TriMesh {
    public:
        GTS_TriMesh( const TriangleMesh<REAL> &mesh );

        // Destructor
        virtual ~GTS_TriMesh();

        // Gaussian curvature at a vertex
        REAL       gaussianCurvature( int vertex_idx );

        // Mesh curvature at a vertex
        Vector3d   meanCurvature( int vertex_idx );

        // Samples mean curvature on a triangle: assumes that mean
        // curvatures have been precomputed
        REAL       meanCurvature( int triangle_idx,
                                  const Vector3d &barycentricPosition ) const;

        // Gets all curvature information at once at a vertex
        void       allCurvatures( int vertex_idx,
                                  REAL &Kg, Vector3d &Kh,
                                  REAL &principalCurvature1,
                                  REAL &principalCurvature2 );

        void precomputeMeanCurvatures( bool smooth = false );

    protected:

    private:
        void initSurface();

        void clear();

        void smoothCurvatures();

    private:
        const TriangleMesh<REAL>    &_mesh;

        // GTS representation of the surface
        GtsSurface                  *_surface;

        vector<GtsVertex *>          _vertices;
        vector<GtsFace *>            _faces;

        // Precomputed mean curvatures for each vertex
        FloatArray                   _meanCurvatures;

};

#endif
