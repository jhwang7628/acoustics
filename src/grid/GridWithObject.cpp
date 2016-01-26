#include "Grid.h"
#include "GridWithObject.h"
#include "signalprocessing/FilterDesign.h" 

void UniformGridWithObject::Reinitialize(const std::string &configFile)
{
    std::cout << "config file : " << configFile << std::endl;

    parser_.reset(Parser::buildParser( configFile )); 
    if ( !parser_ ) throw std::runtime_error("**ERROR** Could not build parser from "+configFile);

    mesh_.reset(parser_->getMesh("impulse_response")); 
    if ( !mesh_ ) throw std::runtime_error("**ERROR** Could not build mesh");

    solverParameters_ = parser_->getImpulseResponseParms();

    distanceField_.reset(DistanceFieldBuilder::BuildSignedClosestPointField(parser_->getMeshFileName().c_str(), solverParameters_._sdfResolution, solverParameters_._sdfFilePrefix.c_str())); 
    if (!distanceField_) throw std::runtime_error("**ERROR** Could not construct distance field"); 

    isBoundary_.clear(); 
    isBoundary_.resize(N_cells());

    interfacialCellTypes_.clear(); 
    interfacialCellTypes_.resize(N_cells()); 

    initialized_ = true; 
}

void UniformGridWithObject::ClassifyCells()
{
    assert(initialized_); 

    std::cout << "classifying cells" << std::endl;
    const Eigen::Vector3i &cellCount = this->cellCount_; 
    const int N_cells = this->N_cells(); 


    // first rasterize the geometry
    int N_solids = 0; 
    for (int kk=0; kk<cellCount[2]; kk++) 
    {
        std::cout << " classify bulk cells. progress: " << kk << "/" << cellCount[2] << "\r" << std::flush;
        for (int jj=0; jj<cellCount[1]; jj++) 
        {
            for (int ii=0; ii<cellCount[0]; ii++) 
            {
                Vector3d position; 
                this->GetCellCenterPosition(ii,jj,kk,position.x,position.y,position.z); 
                if (distanceField_->distance(position) <= distanceTolerance_) 
                {
                    isBoundary_[this->FlattenIndicies(ii,jj,kk)] = SOLID; 
                    N_solids ++; 
                }
            }
        }
    }
    std::cout << std::endl;


    // then find the interfacial cells 
    int index; 
    int N_interfacial = 0; 
    for (int kk=1; kk<cellCount[2]-1; kk++) 
    {
        std::cout << " classify interfacial cells. progress: " << kk << "/" << cellCount[2] << "\r" << std::flush;
        for (int jj=1; jj<cellCount[1]-1; jj++) 
        {
            for (int ii=1; ii<cellCount[0]-1; ii++) 
            {
                index = FlattenIndicies(ii,jj,kk); 
                if (isBoundary_[index]) 
                    continue; 

                InterfacialType &type = interfacialCellTypes_[index]; 

                if (isBoundary_[FlattenIndicies(ii-1,jj  ,kk  )] == SOLID) type = ( type | X_BOUNDARY_ON_LEFT ); 
                if (isBoundary_[FlattenIndicies(ii+1,jj  ,kk  )] == SOLID) type = ( type | X_BOUNDARY_ON_RIGHT);
                if (isBoundary_[FlattenIndicies(ii  ,jj-1,kk  )] == SOLID) type = ( type | Y_BOUNDARY_ON_LEFT ); 
                if (isBoundary_[FlattenIndicies(ii  ,jj+1,kk  )] == SOLID) type = ( type | Y_BOUNDARY_ON_RIGHT); 
                if (isBoundary_[FlattenIndicies(ii  ,jj  ,kk-1)] == SOLID) type = ( type | Z_BOUNDARY_ON_LEFT ); 
                if (isBoundary_[FlattenIndicies(ii  ,jj  ,kk+1)] == SOLID) type = ( type | Z_BOUNDARY_ON_RIGHT); 

                if (type & (X_BOUNDARY_ON_LEFT | X_BOUNDARY_ON_RIGHT | Y_BOUNDARY_ON_LEFT | Y_BOUNDARY_ON_RIGHT | Z_BOUNDARY_ON_LEFT | Z_BOUNDARY_ON_RIGHT) )
                {
                    type = ( type | IS_INTERFACE ); 
                    N_interfacial++; 
                }

            }
        }
    }
    std::cout << std::endl;

    std::cout << "classification completed : \n"; 
    std::cout << " solid : " << N_solids << "\n"; 
    std::cout << " fluid : " << N_cells - N_solids << "\n";
    std::cout << " interfacial cells (neighbor solid) : " << N_interfacial << std::endl;
}

void UniformGridWithObject::WriteCellTypes(const std::string &filename)
{
    std::shared_ptr<Eigen::MatrixXd> cellTypes(new Eigen::MatrixXd(N_cells(),1)); 
    cellTypes->setZero(); 
    Eigen::MatrixXd positionBuffer(N_cells(),3); 
    int indBuffer=0; 
    std::ofstream of(filename.c_str()); 
    for (int kk=0; kk<cellCount_[2]; kk++) 
        for (int jj=0; jj<cellCount_[1]; jj++) 
            for (int ii=0; ii<cellCount_[0]; ii++) 
            {
                FlattenIndicies(ii,jj,kk,indBuffer); 
                GetCellCenterPosition(ii,jj,kk,positionBuffer(indBuffer,0),positionBuffer(indBuffer,1),positionBuffer(indBuffer,2)); 
                if (isBoundary_[indBuffer]) 
                    (*cellTypes)(indBuffer,0) = -1.0;
                if (interfacialCellTypes_[indBuffer] & IS_INTERFACE) 
                    (*cellTypes)(indBuffer,0) =  1.0; 
            }

    InsertCellCenteredData("types", cellTypes); 
    WriteVTKCellCentered("cellTypes", "types", "cellTypes"); 
}


// use first order finite-difference on boundaries 
void UniformGridWithObject::CellCenteredDataGradient( const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &gradientData, const CELL_GRADIENT_COMPONENT &component)
{

    const GridData & data = this->GetCellCenteredData( dataName ); 
    const int dataDimension =  data.NData(); 
    const int NCell = data.NCells(); 
    assert( dataDimension == 1 ); 

    const int Nx = this->cellCount_[0]; 
    const int Ny = this->cellCount_[1]; 
    const int Nz = this->cellCount_[2]; 

    const double dx = ( this->maxBound_[0] - this->minBound_[0] )/(double)Nx; 
    const double dy = ( this->maxBound_[1] - this->minBound_[1] )/(double)Ny; 
    const double dz = ( this->maxBound_[2] - this->minBound_[2] )/(double)Nz; 

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
                double buffer0=0, buffer1=0, buffer2=0; 

                // do the computation only if its non solid
                if (isBoundary_[flattenedInd] != SOLID)
                {
                    int xSign, ySign, zSign; // for boundary 

                    if      (ii==0   )  {xSign =  1;} 
                    else if (ii==Nx-1)  {xSign = -1;} 
                    else                
                    {
                        if      (isBoundary_[FlattenIndicies(ii-1,jj,kk)]==SOLID) xSign =  1;
                        else if (isBoundary_[FlattenIndicies(ii+1,jj,kk)]==SOLID) xSign = -1;
                        else    xSign =  0; // not on boundary 
                    }

                    if      (jj==0)     {ySign =  1;} 
                    else if (jj==Ny-1)  {ySign = -1;} 
                    else                
                    {
                        if      (isBoundary_[FlattenIndicies(ii,jj-1,kk)]==SOLID) ySign =  1;
                        else if (isBoundary_[FlattenIndicies(ii,jj+1,kk)]==SOLID) ySign = -1;
                        else    ySign =  0; // not on boundary 
                    }

                    if      (kk==0)     {zSign =  1;} 
                    else if (kk==Nz-1)  {zSign = -1;} 
                    else               
                    {
                        if      (isBoundary_[FlattenIndicies(ii,jj,kk-1)]==SOLID) zSign =  1;
                        else if (isBoundary_[FlattenIndicies(ii,jj,kk+1)]==SOLID) zSign = -1;
                        else    zSign =  0; // not on boundary 
                    }

                    // second order accurate, but it creates a little jagged
                    // edges since it extended one more stencil to maintain 2nd
                    // accuracy
                    if (xSign != 0) 
                        buffer0 = static_cast<double>(xSign)*(-data.Value(ii,jj,kk)(0) + data.Value(ii+xSign*1,jj,kk)(0))/dx; 
                        //buffer0 = static_cast<double>(xSign)*(-3.*data.Value(ii,jj,kk)(0) + 4.*data.Value(ii+xSign*1,jj,kk)(0) - data.Value(ii+xSign*2,jj,kk)(0))/2./dx; 
                    else 
                        buffer0 = (data.Value(ii+1,jj,kk)(0) - data.Value(ii-1,jj,kk)(0))/2./dx; 

                    if (ySign != 0) 
                        buffer1 = static_cast<double>(ySign)*(-data.Value(ii,jj,kk)(0) + data.Value(ii,jj+ySign*1,kk)(0))/dy; 
                        //buffer1 = static_cast<double>(ySign)*(-3.*data.Value(ii,jj,kk)(0) + 4.*data.Value(ii,jj+ySign*1,kk)(0) - data.Value(ii,jj+ySign*2,kk)(0))/2./dy; 
                    else 
                        buffer1 = (data.Value(ii,jj+1,kk)(0) - data.Value(ii,jj-1,kk)(0))/2./dy; 

                    if (zSign != 0) 
                        buffer2 = static_cast<double>(zSign)*(-data.Value(ii,jj,kk)(0) + data.Value(ii,jj,kk+zSign*1)(0))/dz; 
                        //buffer2 = static_cast<double>(zSign)*(-3.*data.Value(ii,jj,kk)(0) + 4.*data.Value(ii,jj,kk+zSign*1)(0) - data.Value(ii,jj,kk+zSign*2)(0))/2./dz; 
                    else 
                        buffer2 = (data.Value(ii,jj,kk+1)(0) - data.Value(ii,jj,kk-1)(0))/2./dz; 
                }


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

    // local smoothing on interfacial cells to lower the impact on stencil
    // choice on the boundaries
    std::cout << "interface smoothing" << std::endl;
    Eigen::VectorXd filter = SIGNAL_PROCESSING::DiscreteGaussian1D(3,1); 

    std::vector<bool> filterMask(NCell, 1); 
    for (int ii=0; ii<NCell; ii++) 
        filterMask[ii] = ((interfacialCellTypes_[ii]&IS_INTERFACE) ? true : false); 

    for (size_t ii=0; ii<N_gradientNeeded; ii++) 
    {
        InsertCellCenteredData("tmp_g", gradientData[ii]); 
        *(gradientData[ii]) = CellCenteredSmoothing("tmp_g", filter, filterMask); 
        DeleteCellCenteredData("tmp_g"); 
    }




}

