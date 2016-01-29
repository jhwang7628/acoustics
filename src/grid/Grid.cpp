#include "Grid.h" 
#include "vtkConverter/vtkConverter.h"

int GridData::FlattenIndicies(const int x, const int y, const int z) const 
{ 
    return owner_.FlattenIndicies(x,y,z); 
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

void Grid::InsertCellCenteredData( const std::string &name, std::shared_ptr<Eigen::Matrix<REAL,Eigen::Dynamic,Eigen::Dynamic>> dataIn )
{
    if ( dataIn->rows() != cellCount_[0]*cellCount_[1]*cellCount_[2] )
    {
        std::string errmsg  = "**ERROR** wrong cell-centered data size: " + std::to_string(dataIn->rows()) + "x" + std::to_string(dataIn->cols()); 
        errmsg += ", because it is defined on grid of cell count " + std::to_string(cellCount_[0]) + "x"
            + std::to_string(cellCount_[1]) + "x"
            + std::to_string(cellCount_[2]); 
        throw std::runtime_error( errmsg );
    }

    GridData gridData( name, dataIn, *this ); 
    std::pair<std::string,GridData&> mapData (name, gridData);
    cellCenteredData_.insert( mapData );
}

void Grid::InsertVertexData( const std::string &name, std::shared_ptr<Eigen::Matrix<REAL,Eigen::Dynamic,Eigen::Dynamic>> dataIn )
{
    if ( dataIn->rows() != (cellCount_[0])*(cellCount_[1])*(cellCount_[2]) )
    {
        std::string errmsg  = "**ERROR** wrong vertex data size: " + std::to_string(dataIn->rows()) + "x" + std::to_string(dataIn->cols()); 
        errmsg += ", because it is defined on grid of cell count " + std::to_string(cellCount_[0]) + "x"
            + std::to_string(cellCount_[1]) + "x"
            + std::to_string(cellCount_[2]); 
        throw std::runtime_error( errmsg );
    }

    GridData gridData( name, dataIn, *this ); 
    std::pair<std::string,GridData& > mapData (name, gridData);
    vertexData_.insert( mapData );
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////


UniformGrid::UniformGrid(const Eigen::Vector3d& minBound, const Eigen::Vector3d& maxBound, const Eigen::Vector3i& cellCount) : 
    Grid(minBound, maxBound, cellCount)
{
    RecomputeCachedField(); 
}


UniformGrid::UniformGrid(const UniformGrid& grid) : 
    UniformGrid(grid.minBound_,grid.maxBound_,grid.cellCount_) 
{
    RecomputeCachedField();
}

Eigen::Vector3d UniformGrid::GetCellCenterPosition( const int &ii, const int &jj, const int &kk ) const 
{
    Eigen::Vector3d pos; 
    GetCellCenterPosition(ii,jj,kk,pos[0],pos[1],pos[2]);
    return pos;
}
void UniformGrid::GetCellCenterPosition( const int &ii, const int &jj, const int &kk, double &x, double &y, double &z ) const
{
    x = minBound_[0]+((double)(ii)+0.5)*dx_[0]; 
    y = minBound_[1]+((double)(jj)+0.5)*dx_[1]; 
    z = minBound_[2]+((double)(kk)+0.5)*dx_[2]; 
}

void UniformGrid::GetAllCellCenterPosition(Eigen::MatrixXd &gridPosition) const 
{
    gridPosition.resize(N_cells(), 3); 

    int indexBuffer;
    for (int kk=0; kk<cellCount_[2]; kk++) 
        for (int jj=0; jj<cellCount_[1]; jj++) 
            for (int ii=0; ii<cellCount_[0]; ii++) 
            {
                indexBuffer = FlattenIndicies(ii,jj,kk); 
                GetCellCenterPosition(ii,jj,kk,gridPosition(indexBuffer,0),gridPosition(indexBuffer,1),gridPosition(indexBuffer,2));
            }
}

void UniformGrid::GetAllCellCenterPosition(const Eigen::Vector3d &minBound, const Eigen::Vector3d &maxBound, const Eigen::Vector3i &cellCount, Eigen::MatrixXd &gridPosition)
{
    UniformGrid grid(minBound,maxBound,cellCount);
    grid.GetAllCellCenterPosition(gridPosition); 
}

Eigen::VectorXd UniformGrid::InterpolateVertexData( const std::string& dataName, const Eigen::Vector3d& position ) 
{
    const GridData& data = GetVertexData( dataName ); 

    Eigen::VectorXd value; 

    switch (interpolationScheme_)
    {
        case trilinear:

            value = TrilinearInterpolate( minBound_, maxBound_, cellCount_, position, data );

            break; 
    }

    return value; 
}

Eigen::VectorXd UniformGrid::InterpolateCellCenteredData( const std::string& dataName, const Eigen::Vector3d& position )
{

    const GridData& data = GetCellCenteredData( dataName ); 

    Eigen::VectorXd value; 

    switch (interpolationScheme_)
    {
        case trilinear:

            // define a few useful quantities for cell center
            // interpolation
            Eigen::Vector3d cmin = minBound_ + dx_/2.0; 
            Eigen::Vector3d cmax = maxBound_ - dx_/2.0; 
            Eigen::Vector3i Nc = cellCount_ - Eigen::Vector3i::Ones(3); 
            value = TrilinearInterpolate( cmin, cmax, Nc, position, data );

            break; 
    }

    return value; 
}


// this is a naive implementation of the smoothing operation using
// filter by N^3 iteration. 
//
// TODO
// if separable filters are used, the complexity can be brought down by
// using smart caching, see e.g.
// http://blogs.mathworks.com/steve/2006/10/04/separable-convolution/
Eigen::MatrixXd UniformGrid::CellCenteredSmoothing( const std::string &dataName, const Eigen::VectorXd &filter, const std::vector<bool> &mask)
{
    const GridData &data = GetCellCenteredData( dataName ); 
    const int dataDimension = data.NData(); 
    const int N_cells = data.NCells(); 
    assert( dataDimension == 1 ); 
    assert( filter.size() %2 == 1);

    const int kernelHalfWidth = (filter.size()-1)/2;

    const int Nx = cellCount_[0]; 
    const int Ny = cellCount_[1]; 
    const int Nz = cellCount_[2]; 

    Eigen::MatrixXd smoothedData(N_cells,dataDimension); 
    smoothedData.setZero(); 

    int count = 0; 

    int flattenedIndex = 0;

#pragma omp parallel for
    for (int kk=0; kk<Nz; kk++) 
    {
#pragma omp critical
        {
            //std::cout << count << " " << std::flush; 
            count++; 
        }
        for (int jj=0; jj<Ny; jj++) 
        {
            for (int ii=0; ii<Nx; ii++) 
            {
                flattenedIndex = data.FlattenIndicies(ii,jj,kk); 

                if (!mask[flattenedIndex])
                {
                    smoothedData(flattenedIndex,0) = data.Value(ii,jj,kk)(0); 
                    continue; 
                }

                // not smoothing the boundary
                if ((ii-kernelHalfWidth < 0 || ii+kernelHalfWidth>=Nx) || 
                        (jj-kernelHalfWidth < 0 || jj+kernelHalfWidth>=Ny) ||
                        (kk-kernelHalfWidth < 0 || kk+kernelHalfWidth>=Nz))
                {
#pragma omp critical
                    smoothedData(flattenedIndex,0) = data.Value(ii,jj,kk)(0); 
                    continue; 
                }

                double doubleBuffer = 0.0;

                // smooth data
                for (int kk_ker=-kernelHalfWidth; kk_ker<kernelHalfWidth+1; kk_ker++) 
                {
                    for (int jj_ker=-kernelHalfWidth; jj_ker<kernelHalfWidth+1; jj_ker++) 
                    {
                        for (int ii_ker=-kernelHalfWidth; ii_ker<kernelHalfWidth+1; ii_ker++) 
                        {

                            const int kk_ker_shifted = kk_ker+kernelHalfWidth; 
                            const int jj_ker_shifted = jj_ker+kernelHalfWidth; 
                            const int ii_ker_shifted = ii_ker+kernelHalfWidth; 

                            const double filterValue = filter(ii_ker_shifted)*filter(jj_ker_shifted)*filter(kk_ker_shifted); 

                            doubleBuffer += data.Value(ii+ii_ker,jj+jj_ker,kk+kk_ker)(0)*filterValue;
                            //smoothedData(flattenedIndex,0) += data.Value(ii+ii_ker,jj+jj_ker,kk+kk_ker)(0)*filterValue;
                        }
                    }
                }

#pragma omp critical
                smoothedData(flattenedIndex,0) = doubleBuffer;
            }
        }
    }
    //std::cout << std::endl;

    return smoothedData; 

}

void UniformGrid::CellCenteredScalarHessian(const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &hessianData)
{
    std::vector<std::shared_ptr<Eigen::MatrixXd>> gradientBuffer; 
    CellCenteredScalarHessian(dataName, gradientBuffer, hessianData);
}



void UniformGrid::CellCenteredScalarHessian(const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &gradientBuffer, std::vector<std::shared_ptr<Eigen::MatrixXd>> &hessianData)
{

    const GridData &gridData = GetCellCenteredData(dataName); 
    const int dataDimension = gridData.NData(); 
    const int N_cells = gridData.NCells(); 
    assert( dataDimension==1 ); 

    const int hessianComponents = 6;


    //std::vector<std::shared_ptr<Eigen::MatrixXd>> gradientBuffer; 

    CellCenteredDataGradient( dataName, gradientBuffer, ALL); 

    std::shared_ptr<Eigen::MatrixXd> gx = gradientBuffer[0]; 
    std::shared_ptr<Eigen::MatrixXd> gy = gradientBuffer[1]; 
    std::shared_ptr<Eigen::MatrixXd> gz = gradientBuffer[2]; 

    InsertCellCenteredData(dataName+"gx", gx);
    InsertCellCenteredData(dataName+"gy", gy);
    InsertCellCenteredData(dataName+"gz", gz);

    std::vector<std::shared_ptr<Eigen::MatrixXd>> d_gx; 
    std::vector<std::shared_ptr<Eigen::MatrixXd>> d_gy; 
    std::vector<std::shared_ptr<Eigen::MatrixXd>> d_gz; 

    d_gx.push_back(hessianData[0]); 
    d_gx.push_back(hessianData[1]); 
    d_gx.push_back(hessianData[2]); 
    d_gy.push_back(hessianData[3]); 
    d_gy.push_back(hessianData[4]); 
    d_gz.push_back(hessianData[5]); 

    CellCenteredDataGradient(dataName+"gx", d_gx, ALL); 
    CellCenteredDataGradient(dataName+"gy", d_gy, YZ ); 
    CellCenteredDataGradient(dataName+"gz", d_gz, Z  ); 

    bool resizeNeeded = false;
    if (hessianData.size()!=hessianComponents)
        resizeNeeded = true; 
    else 
    {
        for (auto &p : hessianData) 
        {
            if (p == nullptr) { resizeNeeded = true; } 
            else { if (p->rows()!=N_cells || p->cols()!=dataDimension) {resizeNeeded = true; } }
        }
    }

    if (resizeNeeded) 
    {
        throw std::runtime_error("**ERROR** should not resize in hessian");
        hessianData.clear(); 
        hessianData.resize(6); 
    }

    DeleteCellCenteredData(dataName+"gx"); 
    DeleteCellCenteredData(dataName+"gy"); 
    DeleteCellCenteredData(dataName+"gz"); 

}


// 2nd finite-difference for 1st-order derivative 
//   (f_i+1 - f_i-1) / 2h,              if 0<i<N-1
//   (-3f_i + 4f_i+1 - f_i+2) / 2h,     if i=0
//   (+3f_i - 4f_i-1 + f_i-2) / 2h,     if i=N-1
//
void UniformGrid::CellCenteredDataGradient( const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &gradientData , const CELL_GRADIENT_COMPONENT &component)
{
    const GridData & data = GetCellCenteredData( dataName ); 
    const int dataDimension =  data.NData(); 
    const int NCell = data.NCells(); 
    assert( dataDimension == 1 ); 

    const int Nx = cellCount_[0]; 
    const int Ny = cellCount_[1]; 
    const int Nz = cellCount_[2]; 

    const double dx = ( maxBound_[0] - minBound_[0] )/(double)Nx; 
    const double dy = ( maxBound_[1] - minBound_[1] )/(double)Ny; 
    const double dz = ( maxBound_[2] - minBound_[2] )/(double)Nz; 

    size_t N_gradientNeeded; 
    const int i_component = static_cast<int>(component); 
    if (i_component==0)
        N_gradientNeeded = 3; 
    else if (i_component>0 && i_component<10)
        N_gradientNeeded = 1; 
    else if (i_component>10)
        N_gradientNeeded = 2; 
    else 
        throw std::runtime_error("**ERROR** wrong cell ponent input : "+std::to_string(i_component)); 

    // determine whether to resize the input vector
    bool resizeNeeded = false; 
    for (auto & m : gradientData) 
    {
        if (nullptr == m) 
        {
            resizeNeeded = true; 
            break; 
        }

        if (m->rows()!=NCell || m->cols()!=1)
        {
            resizeNeeded = true; 
            break; 
        }
    }
    if (gradientData.size()!=N_gradientNeeded)
        resizeNeeded = true; 

    if (resizeNeeded) 
    {
        throw std::runtime_error("**ERROR** should not resize in gradient");
        gradientData.clear(); 
        gradientData.resize(N_gradientNeeded); 
        for (size_t ii=0; ii<N_gradientNeeded; ii++) 
            gradientData[ii].reset(new Eigen::MatrixXd(NCell,1));
    }

    // loop through each direction and compute the data difference and
    // take care of the boundaries
#pragma omp parallel for
    for (int kk=0; kk<Nz; kk++)
    {
        for (int jj=0; jj<Ny; jj++)
        {
            for (int ii=0; ii<Ny; ii++)
            { 
                const int flattenedInd = data.FlattenIndicies(ii,jj,kk); 

                int xSign, ySign, zSign; // for boundary 

                if      (ii==0   )  xSign =  1; 
                else if (ii==Nx-1)  xSign = -1; 
                else                xSign =  0; // not on boundary 

                if      (jj==0)     ySign =  1; 
                else if (jj==Ny-1)  ySign = -1; 
                else                ySign =  0; // not on boundary 

                if      (kk==0)     zSign =  1; 
                else if (kk==Nz-1)  zSign = -1; 
                else                zSign =  0; // not on boundary 



                double buffer0, buffer1, buffer2; 
                if (xSign != 0) 
                    buffer0 = static_cast<double>(xSign)*(-3.*data.Value(ii,jj,kk)(0) + 4.*data.Value(ii+xSign*1,jj,kk)(0) - data.Value(ii+xSign*2,jj,kk)(0))/2./dx; 
                else 
                    buffer0 = (data.Value(ii+1,jj,kk)(0) - data.Value(ii-1,jj,kk)(0))/2./dx; 

                if (ySign != 0) 
                    buffer1 = static_cast<double>(ySign)*(-3.*data.Value(ii,jj,kk)(0) + 4.*data.Value(ii,jj+ySign*1,kk)(0) - data.Value(ii,jj+ySign*2,kk)(0))/2./dy; 
                else 
                    buffer1 = (data.Value(ii,jj+1,kk)(0) - data.Value(ii,jj-1,kk)(0))/2./dy; 

                if (zSign != 0) 
                    buffer2 = static_cast<double>(zSign)*(-3.*data.Value(ii,jj,kk)(0) + 4.*data.Value(ii,jj,kk+zSign*1)(0) - data.Value(ii,jj,kk+zSign*2)(0))/2./dz; 
                else 
                    buffer2 = (data.Value(ii,jj,kk+1)(0) - data.Value(ii,jj,kk-1)(0))/2./dz; 


#pragma omp critical
                {
                    switch (i_component)
                    {
                        case 0: 
                            (*(gradientData[0]))(flattenedInd, 0) = buffer0; 
                            (*(gradientData[1]))(flattenedInd, 0) = buffer1; 
                            (*(gradientData[2]))(flattenedInd, 0) = buffer2; 
                            break;
                        case 1: 
                            (*(gradientData[0]))(flattenedInd, 0) = buffer0; 
                            break;
                        case 2: 
                            (*(gradientData[0]))(flattenedInd, 0) = buffer1; 
                            break;
                        case 3: 
                            (*(gradientData[0]))(flattenedInd, 0) = buffer2; 
                            break;
                        case 12: 
                            (*(gradientData[0]))(flattenedInd, 0) = buffer0; 
                            (*(gradientData[1]))(flattenedInd, 0) = buffer1; 
                            break;
                        case 23: 
                            (*(gradientData[0]))(flattenedInd, 0) = buffer1; 
                            (*(gradientData[1]))(flattenedInd, 0) = buffer2; 
                            break;
                        case 13: 
                            (*(gradientData[0]))(flattenedInd, 0) = buffer0; 
                            (*(gradientData[1]))(flattenedInd, 0) = buffer2; 
                            break;
                    }
                }
            } 
        } 
    }

}

void UniformGrid::PrintDataMeta()
{
    std::cout << "Print all data on the grid:\n"; 
    std::cout << "---------------------------\n"; 


    std::cout << "Cell-centered:\n";
    std::cout << cellCenteredData_.size() << std::endl;
}

void UniformGrid::WriteVTKCellGrid(const std::string &vtkName)
{
    const Eigen::Vector3i cellCount = GetCellCount(); 
    Eigen::MatrixXd gridPosition(cellCount[0]*cellCount[1]*cellCount[2],3); 

    for (int kk=0; kk<cellCount[2]+1; kk++) 
    {
        for (int jj=0; jj<cellCount[1]+1; jj++)
        {
            for (int ii=0; ii<cellCount[0]+1; ii++)
            {
                const int index  = FlattenIndicies(ii,jj,kk); 
                gridPosition.row(index) = GetCellCenterPosition(ii,jj,kk)-dx_/2.;
            }
        }
    }

    VTKConverter::VTKStructureGridFromEigen( gridPosition, vtkName, VTKConverter::BINARY, cellCount );
}

void UniformGrid::WriteVTKCellCentered(const std::string &vtkPrefix, const std::string &key, const std::string &dataName)
{
    const Eigen::Vector3i cellCount = GetCellCount(); 
    Eigen::MatrixXd gridPosition(cellCount[0]*cellCount[1]*cellCount[2],3); 
    GridData &gridData = GetCellCenteredData(key);
    Eigen::MatrixXd &data = gridData.GetMatrix(); 

    const int N_data = data.cols(); 

    for (int kk=0; kk<cellCount[2]; kk++) 
    {
        for (int jj=0; jj<cellCount[1]; jj++)
        {
            for (int ii=0; ii<cellCount[0]; ii++)
            {
                const int index  = gridData.FlattenIndicies(ii,jj,kk); 
                gridPosition.row(index) = GetCellCenterPosition(ii,jj,kk);
            }
        }
    }

    if (N_data>1)
    {
        for (int ii=0; ii<N_data; ii++) 
            VTKConverter::VTKStructureGridWithScalarFromEigen( gridPosition, data.col(ii), vtkPrefix+"_"+std::to_string(ii)+".vtk", dataName, VTKConverter::BINARY, cellCount );
    }
    else 
        VTKConverter::VTKStructureGridWithScalarFromEigen( gridPosition, data.col(0), vtkPrefix+".vtk", dataName, VTKConverter::BINARY, cellCount );
}

void UniformGrid::WriteVTKCellCenteredFromEigen(std::shared_ptr<Eigen::MatrixXd> ptrData, std::shared_ptr<Eigen::MatrixXd> ptrGridPosition, const std::string &vtkPrefix, const std::string &dataName)
{
}



void UniformGrid::WriteVTKCellCenteredFromEigen(std::shared_ptr<Eigen::MatrixXd> ptrData, const Eigen::Vector3d &minBound, const Eigen::Vector3d &maxBound, const Eigen::Vector3i &cellCount, const std::string &vtkPrefix, const std::string &dataName)
{

    UniformGrid grid(minBound,maxBound,cellCount); 
    const std::string tmpKey("tmp_data"); 
    grid.InsertCellCenteredData(tmpKey, ptrData); 
    grid.WriteVTKCellCentered(vtkPrefix, tmpKey, dataName);
    grid.DeleteCellCenteredData(tmpKey); 
}

void UniformGrid::WriteVTKCellCenteredFromEigen(const Eigen::MatrixXd &data, const Eigen::Vector3d &minBound, const Eigen::Vector3d &maxBound, const Eigen::Vector3i &cellCount, const std::string &vtkPrefix, const std::string &dataName)
{
    std::shared_ptr<Eigen::MatrixXd> ptrData(new Eigen::MatrixXd(data)); 
    UniformGrid::WriteVTKCellCenteredFromEigen(ptrData, minBound, maxBound, cellCount, vtkPrefix, dataName); 
}

void UniformGrid::WriteCSVCellCentered(const std::string &csvName, const std::string &dataName)
{
    std::ofstream of(csvName.c_str());
    const int &Nx = cellCount_[0]; 
    const int &Ny = cellCount_[1]; 
    const int &Nz = cellCount_[2]; 
    const Eigen::Vector3d &maxBound = maxBound_;
    const Eigen::Vector3d &minBound = minBound_;
    const Eigen::Vector3i &division = cellCount_;
    Eigen::MatrixXd &data = GetCellCenteredData(dataName).GetMatrix(); 

    for (int ii=0; ii<Nx; ii++) 
    {
        for (int jj=0; jj<Ny; jj++) 
        {
            for (int kk=0; kk<Nz; kk++) 
            {
                const int ind = ii + jj*Nx + kk*Nx*Ny; 
                const double x = minBound[0] + ii*(maxBound[0]-minBound[0])/division[0] + dx_[0]/2.; 
                const double y = minBound[1] + jj*(maxBound[1]-minBound[1])/division[1] + dx_[1]/2.; 
                const double z = minBound[2] + kk*(maxBound[2]-minBound[2])/division[2] + dx_[2]/2.; 
                of << x << "," << y << "," << z << "," << data(ind,0) << std::endl;
            }
        }
    }
    of.close();
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
void RotatedUniformGrid::SetRotationMatrix()
{
    Eigen::Matrix3d rotation; 
    rotation << cos(yaw_), 0,       sin(yaw_), 
             0        , 1,       0, 
             -sin(yaw_), 0,       cos(yaw_); 
    rotation_ = rotation; 
}

Eigen::Vector3d RotatedUniformGrid::GetCellCenterPosition( const int ii, const int jj, const int kk ) const 
{

    Eigen::Vector3d position(minBound_[0]+((double)(ii)+0.5)*dx_[0], 
            minBound_[1]+((double)(jj)+0.5)*dx_[1],
            minBound_[2]+((double)(kk)+0.5)*dx_[2]); 

    return rotation_*position; 

}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
/*
 * Helper for trilinear interpolation on the grid.
 *
 * trilinear interpolation following notation on wiki
 * (https://en.wikipedia.org/wiki/Trilinear_interpolation)
 *
 * this implementation does not depend on grid size
 */
Eigen::VectorXd TrilinearInterpolate( const Eigen::Vector3d& xmin, const Eigen::Vector3d& xmax, const Eigen::Vector3i& N, const Eigen::Vector3d& position, const GridData& data ) 
{
    Eigen::VectorXd value; 
    value.setZero(data.NData());

    // get bounding index
    Eigen::Vector3i lowIndex; 
    Eigen::Vector3i highIndex;
    Eigen::Vector3d dx = (xmax - xmin).cwiseQuotient(N.cast<double>());
    for (int ii=0; ii<3; ii++) 
    {
        lowIndex[ii] = std::min<int>( floor((position[ii]-xmin[ii])/dx[ii]), N[ii]-1 );
        lowIndex[ii] = std::max<int>( lowIndex[ii], 0 ); // prevent rounding error
        highIndex[ii] = lowIndex[ii] + 1; 
    }

    int ix0,iy0,iz0,ix1,iy1,iz1; 
    ix0= lowIndex[0]; 
    ix1=highIndex[0]; 
    iy0= lowIndex[1]; 
    iy1=highIndex[1]; 
    iz0= lowIndex[2]; 
    iz1=highIndex[2]; 

    double x0,y0,z0,x1,y1,z1,xd,yd,zd,x,y,z; 
    x = std::max<REAL>( std::min<REAL>( position[0], xmax[0] ), xmin[0] );
    y = std::max<REAL>( std::min<REAL>( position[1], xmax[1] ), xmin[1] );
    z = std::max<REAL>( std::min<REAL>( position[2], xmax[2] ), xmin[2] );
    x0 = xmin[0] + ix0*dx[0]; 
    y0 = xmin[1] + iy0*dx[1]; 
    z0 = xmin[2] + iz0*dx[2]; 
    x1 = x0 + dx[0]; 
    y1 = y0 + dx[1]; 
    z1 = z0 + dx[2]; 
    xd = (x-x0)/(x1-x0); 
    yd = (y-y0)/(y1-y0); 
    zd = (z-z0)/(z1-z0); 

    for (int ii=0; ii<data.NData(); ii++) 
    {

        REAL c00 = data.Value(ix0,iy0,iz0)(ii)*(1.0-xd) + data.Value(ix1,iy0,iz0)(ii)*xd; 
        REAL c10 = data.Value(ix0,iy1,iz0)(ii)*(1.0-xd) + data.Value(ix1,iy1,iz0)(ii)*xd; 
        REAL c01 = data.Value(ix0,iy0,iz1)(ii)*(1.0-xd) + data.Value(ix1,iy0,iz1)(ii)*xd; 
        REAL c11 = data.Value(ix0,iy1,iz1)(ii)*(1.0-xd) + data.Value(ix1,iy1,iz1)(ii)*xd; 

        REAL c0 = c00*(1.0 - yd) + c10*yd; 
        REAL c1 = c01*(1.0 - yd) + c11*yd; 

        value[ii] = c0*(1.0 - zd) + c1*zd; 
    }

    return value;
}
