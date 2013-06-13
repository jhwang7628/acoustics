#ifndef STVK_MESH_H
#   define STVK_MESH_H

#include "stvk.h"
#include "geometry/FixVtxTetMesh.hpp"
#include "linearalgebra/PardisoMatrix.hpp"

class StVKMesh : public StVKMaterial
{
    public:
        typedef PardisoMatrix<double>   TStiffMat;

        StVKMesh(FixVtxTetMesh<double>& msh, double lam, double mu):
                StVKMaterial(lam, mu), refmesh_(msh)
        { 
            generate_sparse_pattern(stiffMat_); 
        }

        void generate_sparse_pattern(TStiffMat& mat) const;
        const TStiffMat& update_stiffness_matrix();
        const TStiffMat& mesh_stiffness_matrix() const
        {   return stiffMat_; }

    private:
        FixVtxTetMesh<double>&      refmesh_;

        TStiffMat                   stiffMat_;
};

#endif
