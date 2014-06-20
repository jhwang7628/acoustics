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
#include <QMessageBox>
#include <QKeyEvent>
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
    if ( !mesh_ ) {
        // Can't load modes if we don't have a mesh
        QMessageBox msgBox;
        msgBox.setText("Must load a tet mesh first");
        msgBox.setIcon(QMessageBox::Information);
        msgBox.exec();

        return;
    }

    QString file = QFileDialog::getOpenFileName(this,
            "Select the mode file", ".", 
            "mode data binary file (*.modes)");
    if ( file.isEmpty() ) return;

    modeData_.read( file.toAscii().data() );

    if ( modeData_.numDOF() != (int)vtx_.size() * 3 ) {
        QMessageBox msgBox;
        msgBox.setText("Error: Mode size does not match mesh size");
        msgBox.setIcon(QMessageBox::Information);
        msgBox.exec();

        return;
    }

    modeIndex->setEnabled( true );
    modalCoordinate->setEnabled( true );
    modeScale->setEnabled( true );

    modeIndex->setMinimum( 1 );
    modeIndex->setMaximum( modeData_.numModes() );
    modeIndex->setValue( 1 );

    modalCoordinate->setValue( 0 );

    modeScale->setValue( 1.0 );

    activeMode_ = 0;
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
    for(size_t i = 0;i < tgl.size();++ i)
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
    for(size_t i = 0;i < tgl.size();++ i)
    {
        const Vector3d& n = tglnmls[i];
        nml_[tgl[i][0]] += n * Triangle<double>::angle(
                vtx_[tgl[i][2]], vtx_[tgl[i][0]], vtx_[tgl[i][1]]);
        nml_[tgl[i][1]] += n * Triangle<double>::angle(
                vtx_[tgl[i][0]], vtx_[tgl[i][1]], vtx_[tgl[i][2]]);
        nml_[tgl[i][2]] += n * Triangle<double>::angle(
                vtx_[tgl[i][1]], vtx_[tgl[i][2]], vtx_[tgl[i][0]]);
    }

    for(size_t i = 0;i < nml_.size();++ i) {
        if ( nml_[i].lengthSqr() > 1E-14 ) nml_[i].normalize();
    }
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
    for(size_t i = 0;i < idx.size();++ i)
    {
        used[idx[i][0]] = true;
        used[idx[i][1]] = true;
        used[idx[i][2]] = true;
        used[idx[i][3]] = true;
    }

    for(size_t i = 0;i < mesh_->num_vertices();++ i) {
        if ( !used[i] ) {
            PRINT_WARNING("Free vertex that is not in use is detected! VID=%d\n", i);
            return;
        }
    }
    PRINT_MSG("No Free Vertex detected\n");
}

void TetViewerFrame::update_active_mode()
{
    modalCoordinate->setValue( 0 );

    int idx = modeIndex->value() - 1;
    REAL frequency = sqrt( modeData_.omegaSquared( idx ) / objectDensity->value() );
    frequency /= 2.0 * M_PI;
    canvas->update_mode_info( modeIndex->value(), frequency );

    update_mode_displacement();
}

void TetViewerFrame::update_mode_displacement()
{
    const vector<Point3d>& vs = mesh_->vertices();
    memcpy(&vtx_[0], &vs[0], sizeof(Point3d)*vs.size());

    REAL displacementScale = (REAL)modalCoordinate->value() / 100.0;
    displacementScale *= (REAL)modeScale->value();
    const vector<REAL> &mode = modeData_.mode( modeIndex->value() - 1 );

    REAL *vertexData = (REAL *)&vtx_[0];
    printf("VS: %d Mode: %d\n", vs.size(), mode.size());
    for ( int i = 0; i < mode.size(); i++ ) {
        vertexData[i] += displacementScale * mode[i];
    }
    update_normals();

    canvas->updateGL();
}

// ============================================================================
int main(int argc, char* argv[])
{
    QApplication application(argc, argv);

    TetViewerFrame mainwnd;
    mainwnd.show();

    return application.exec();
}
