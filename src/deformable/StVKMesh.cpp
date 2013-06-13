#include "StVKMesh.h"
#include "linearalgebra/Matrix.hpp"

void StVKMesh::generate_sparse_pattern(TStiffMat& mat) const
{
    using namespace std;

    const int NUM_FIXED = refmesh_.num_fixed_vertices();
    const int NUM_FREE  = refmesh_.num_free_vertices();

    mat.resize(NUM_FREE*3, true);
    const vector< Tet<double> >& tets = refmesh_.tets();
    const vector<TetMesh<double>::TetIdx>& idx = refmesh_.tet_indices();
    for(size_t tetid = 0;tetid < tets.size();++ tetid)
    for(size_t nid = 0;nid < 4;++ nid)
    {
        if ( refmesh_.is_fixed_vertex(idx[tetid][nid]) ) continue;

        int rowId = (idx[tetid][nid] - NUM_FIXED)*3;
        for(size_t n2id = 0;n2id < 4;++ n2id)
        {
            if ( refmesh_.is_fixed_vertex(idx[tetid][n2id]) ) continue;

            int colId = (idx[tetid][n2id] - NUM_FIXED)*3;
            for(size_t dir = 0;dir < 3;++ dir)
            for(size_t dir2 = 0;dir2 < 3;++ dir2)
                mat.set_nonzero(rowId+dir, colId+dir2);
        }
    }
    mat.generate_pattern();
}

const StVKMesh::TStiffMat& StVKMesh::update_stiffness_matrix()
{
    using namespace std;

    Matrix<double> stiff(12, true);
    const int NUM_FIXED = refmesh_.num_fixed_vertices();

    stiffMat_.zeros();
    const vector< Tet<double> >& tets = refmesh_.tets();
    const vector<TetMesh<double>::TetIdx>& idx = refmesh_.tet_indices();
    for(size_t tetid = 0;tetid < tets.size();++ tetid)
    {
        stiffness_matrix(tets[tetid], stiff);
        for(size_t nid = 0;nid < 4;++ nid)          // iterate on each node in the tet
        {
            if ( refmesh_.is_fixed_vertex(idx[tetid][nid]) ) continue;

            int rowId = (idx[tetid][nid] - NUM_FIXED)*3;
            for(size_t n2id = 0;n2id < 4;++ n2id)   // each node in the same tet
            {
                if ( refmesh_.is_fixed_vertex(idx[tetid][n2id]) ||
                     idx[tetid][nid] > idx[tetid][n2id] ) continue;

                int colId = (idx[tetid][n2id] - NUM_FIXED)*3;
                for(size_t dir = 0;dir < 3;++ dir)
                for(size_t dir2 = 0;dir2 < 3;++ dir2)
                {
                    int srow = nid*3 + dir;
                    int scol = n2id*3 + dir2;
                    stiffMat_.add(rowId+dir, colId+dir2,
                            srow <= scol ? stiff[srow][scol] : stiff[scol][srow]);
                }
            }
        }
    }

    return stiffMat_;
}

