/*
 * =====================================================================================
 *
 *       Filename:  TetViewerFrame.cpp
 *
 *        Version:  1.0
 *        Created:  11/17/10 16:43:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#include "TetViewerFrame.h"
#include <QFileDialog>
#include "io/StellarIO.hpp"
#include "io/TetMeshReader.hpp"
#include "io/TetMeshWriter.hpp"
#include "utils/print_msg.h"

using namespace std;

void TetViewerFrame::open()
{
    QString file = QFileDialog::getOpenFileName(this,
            "Select the mesh file", ".", 
            "mesh bin file (*.tet);;mesh txt file (*.node)");
    if ( file.isEmpty() ) return;
    
    TMesh* msh = new TMesh;
    if ( file.endsWith(".tet", Qt::CaseInsensitive) )
        FV_TetMeshLoader_Double::load_mesh(file.toAscii().data(), *msh);
        //TetMeshLoader_Double::load_mesh(file.toAscii().data(), *msh);
    else if ( file.endsWith(".node", Qt::CaseInsensitive) )
    {
        StellarTetMeshLoader::load_mesh(file.left(file.length()-5).toAscii().data(), msh);
    }
    msh->init();
    msh->update_surface();

    update_mesh(msh);
}

void TetViewerFrame::load_modes()
{
    QString file = QFileDialog::getOpenFileName(this,
            "Select the mode file", ".", 
            "mode data binary file (*.modes)");
    if ( file.isEmpty() ) return;
    
    TMesh* msh = new TMesh;
    if ( file.endsWith(".tet", Qt::CaseInsensitive) )
        FV_TetMeshLoader_Double::load_mesh(file.toAscii().data(), *msh);
        //TetMeshLoader_Double::load_mesh(file.toAscii().data(), *msh);
    else if ( file.endsWith(".node", Qt::CaseInsensitive) )
    {
        StellarTetMeshLoader::load_mesh(file.left(file.length()-5).toAscii().data(), msh);
    }
    msh->init();
    msh->update_surface();

    update_mesh(msh);
}

void TetViewerFrame::export_bin_tet()
{
    if ( !mesh_ ) return;

    QString file = QFileDialog::getSaveFileName(this,
            "Binary Mesh file name", ".", "bin tet file (*.tet);;All (*)");
    if ( file.isEmpty() ) return;
    FV_TetMeshWriter_Double::write_mesh(file.toAscii().data(), *mesh_);
    //TetMeshWriter_Double::write_mesh(file.toAscii().data(), *mesh_);
}

void TetViewerFrame::export_abaqus_tet()
{
    if ( !mesh_ ) return;
    QString file = QFileDialog::getSaveFileName(this,
            "Abaqus Mesh file name", ".", "abaqus tet file (*.aba);;All (*)");
    if ( file.isEmpty() ) return;
    FV_AbaqusMeshWriter::write_mesh(file.toAscii().data(), *mesh_);
}

void TetViewerFrame::update_mesh(TMesh* msh)
{
    if ( mesh_ )
    {
        delete mesh_;
    }

    vtx_.resize(msh->num_vertices());
    const vector<Point3d>& vs = msh->vertices();
    memcpy(&vtx_[0], &vs[0], sizeof(Point3d)*vs.size());
    mesh_ = msh;

    PRINT_MSG("Update normals\n");
    update_normals();
}

void TetViewerFrame::update_normals()
{
    nml_.resize(mesh_->num_vertices());
    const vector<Tuple3ui>& tgl = mesh_->surface_indices();
    vector<Vector3d> tglnmls(tgl.size());
    for(int i = 0;i < tgl.size();++ i)
    {
        tglnmls[i] = Triangle<double>::normal(
                vtx_[tgl[i][0]], vtx_[tgl[i][1]], vtx_[tgl[i][2]]);
        if ( tglnmls[i].lengthSqr() < 1E-24 )
        {
            PRINT_ERROR("triangle has zero area: %.30lf\n",
                    tglnmls[i].lengthSqr());
            exit(1);
        }
        tglnmls[i].normalize();
    }

    memset(&nml_[0], 0, sizeof(double)*nml_.size());
    for(int i = 0;i < tgl.size();++ i)
    {
        const Vector3d& n = tglnmls[i];
        nml_[tgl[i][0]] += n * Triangle<double>::angle(
                vtx_[tgl[i][2]], vtx_[tgl[i][0]], vtx_[tgl[i][1]]);
        nml_[tgl[i][1]] += n * Triangle<double>::angle(
                vtx_[tgl[i][0]], vtx_[tgl[i][1]], vtx_[tgl[i][2]]);
        nml_[tgl[i][2]] += n * Triangle<double>::angle(
                vtx_[tgl[i][1]], vtx_[tgl[i][2]], vtx_[tgl[i][0]]);
    }

    for(size_t i = 0;i < nml_.size();++ i) 
        if ( nml_[i].lengthSqr() > 1E-14 ) nml_[i].normalize();
}

void TetViewerFrame::check_useless_vtx()
{
    if ( !mesh_ ) 
    {
        PRINT_WARNING("No mesh is loaded yet\n");
        return;
    }

    bool * used = new bool[ mesh_->num_vertices() ];
    memset(used, false, sizeof(bool)*mesh_->num_vertices());

    const vector<TetMesh<double>::TetIdx>& idx = mesh_->tet_indices();
    for(int i = 0;i < idx.size();++ i)
    {
        used[idx[i][0]] = true;
        used[idx[i][1]] = true;
        used[idx[i][2]] = true;
        used[idx[i][3]] = true;
    }

    for(int i = 0;i < mesh_->num_vertices();++ i)
        if ( !used[i] ) 
        {
            PRINT_WARNING("Free vertex that is not in use is detected! VID=%d\n", i);
            return;
        }
    PRINT_MSG("No Free Vertex detected\n");
}

// ============================================================================
int main(int argc, char* argv[])
{
    QApplication application(argc, argv);

    TetViewerFrame mainwnd;
    mainwnd.show();

    return application.exec();
}
