#include "Grid.h"
#include "GridWithObject.h"
#include "signalprocessing/FilterDesign.h" 
#include "geometry/AffineTransformation.h" 
#include "geometry/BoundingBox.h" 

#include "IO/IO.h"


// recover the geometry parameters from the wave solver config file 
void UniformGridWithObject::Reinitialize(const std::string &configFile)
{
    std::cout << "config file : " << configFile << std::endl;

    // parse config file to get geometry parameters 
    parser_.reset(Parser::buildParser( configFile )); 
    if ( !parser_ ) throw std::runtime_error("**ERROR** Could not build parser from "+configFile);

    mesh_.reset(parser_->getMesh("impulse_response")); 
    if ( !mesh_ ) throw std::runtime_error("**ERROR** Could not build mesh");

    solverParameters_ = parser_->getImpulseResponseParms();

    distanceField_.reset(DistanceFieldBuilder::BuildSignedClosestPointField(parser_->getMeshFileName().c_str(), solverParameters_._sdfResolution, solverParameters_._sdfFilePrefix.c_str())); 
    if (!distanceField_) throw std::runtime_error("**ERROR** Could not construct distance field"); 

    cellCount_[0] = solverParameters_._gridResolution; 
    cellCount_[1] = solverParameters_._gridResolution; 
    cellCount_[2] = solverParameters_._gridResolution; 

    BoundingBox fieldBBox(distanceField_->bmin(), distanceField_->bmax()); 
    fieldBBox *= solverParameters_._gridScale; 

    double cellSize = fieldBBox.minlength()/(REAL)solverParameters_._gridResolution; 

    minBound_ = Eigen::Vector3d(fieldBBox.minBound().x, fieldBBox.minBound().y, fieldBBox.minBound().z);
    maxBound_ = minBound_ + Eigen::Vector3d::Ones()*cellSize*(double)solverParameters_._gridResolution; 

    RecomputeCachedField(); 

    // reinitialize fields 

    cellTypes_.clear(); 
    cellTypes_.resize(N_cells(),0); 


    // output 

    std::cout << *this << std::endl;

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
    Vector3d position; 
    for (int kk=0; kk<cellCount[2]; kk++) 
    {
        std::cout << " classify bulk cells. progress: " << kk << "/" << cellCount[2]-1 << "\r" << std::flush;
        for (int jj=0; jj<cellCount[1]; jj++) 
        {
            for (int ii=0; ii<cellCount[0]; ii++) 
            {
                GetCellCenterPosition(ii,jj,kk,position.x,position.y,position.z); 
                if (distanceField_->distance(position) <= distanceTolerance_) 
                {
                    cellTypes_[FlattenIndicies(ii,jj,kk)] |= IS_SOLID; 
                    //isBoundary_[this->FlattenIndicies(ii,jj,kk)] = SOLID; 
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
        std::cout << " classify interfacial cells. progress: " << kk << "/" << cellCount[2]-2 << "\r" << std::flush;
        for (int jj=1; jj<cellCount[1]-1; jj++) 
        {
            for (int ii=1; ii<cellCount[0]-1; ii++) 
            {
                index = FlattenIndicies(ii,jj,kk); 

                if ((cellTypes_[index]&IS_SOLID) == 0) // skip if not solid
                    continue; 

                if (cellTypes_[FlattenIndicies(ii-1,jj  ,kk  )] & IS_SOLID) cellTypes_[index] |= X_SOLID_ON_LEFT ; 
                if (cellTypes_[FlattenIndicies(ii+1,jj  ,kk  )] & IS_SOLID) cellTypes_[index] |= X_SOLID_ON_RIGHT;
                if (cellTypes_[FlattenIndicies(ii  ,jj-1,kk  )] & IS_SOLID) cellTypes_[index] |= Y_SOLID_ON_LEFT ; 
                if (cellTypes_[FlattenIndicies(ii  ,jj+1,kk  )] & IS_SOLID) cellTypes_[index] |= Y_SOLID_ON_RIGHT; 
                if (cellTypes_[FlattenIndicies(ii  ,jj  ,kk-1)] & IS_SOLID) cellTypes_[index] |= Z_SOLID_ON_LEFT ; 
                if (cellTypes_[FlattenIndicies(ii  ,jj  ,kk+1)] & IS_SOLID) cellTypes_[index] |= Z_SOLID_ON_RIGHT; 

                if ((cellTypes_[index] & X_SOLID_ON_LEFT)==0 || (cellTypes_[index] & X_SOLID_ON_RIGHT)==0 || 
                        (cellTypes_[index] & Y_SOLID_ON_LEFT)==0 || (cellTypes_[index] & Y_SOLID_ON_RIGHT)==0 || 
                        (cellTypes_[index] & Z_SOLID_ON_LEFT)==0 || (cellTypes_[index] & Z_SOLID_ON_RIGHT)==0)
                {
                    cellTypes_[index] |= IS_INTERFACE; 
                    N_interfacial++; 
                }

            }
        }
    }
    std::cout << std::endl;

    std::cout << "classification completed : \n"; 
    std::cout << " solid : " << N_solids << "\n"; 
    std::cout << " fluid : " << N_cells - N_solids << "\n";
    std::cout << " interfacial cells (subset of solid cells neighboring fluid) : " << N_interfacial << std::endl;
}

// Find the nearest cell index by even extension off the object boundary,
// described by the sdf and normal. 
//
// it will return the cell index directly if sdf is evaluated to be positive. 
//
// in the negative case, the routine will take normal and closest point, both
// returned by the sdf class, and do certain number of reflections to try to
// get this point out of the object. Once outside, it will then find a cell
// index that is closest to that point which is not classified as solid
// and return the flattened index of that cell. this is essentially a piecewise 
// constant approximation for the zero neumann bc. 
//
// TODO 
// the reflection gives reasonably well approximation
// some problems were found during testing for concave objects. 
//
// 1. sdf uses float internally and the precision could cause rounding error
//    and points being misclassified 
// 2. there is a doubt about how well is the reflection based on the closest
//    point on the geometry and the normal computed using trilinear
//    interpolation. it might be the case that this error converges sublinearly
//    and causes misclassification problem described in 1. 
bool UniformGridWithObject::FlattenIndiciesWithReflection(const int &ii, const int &jj, const int &kk, Eigen::Vector3d &reflectedPositionOut)
{

    int indexBuffer; 
    FlattenIndicies(ii,jj,kk,indexBuffer); 
    GetCellCenterPosition(ii,jj,kk,reflectedPositionOut[0],reflectedPositionOut[1],reflectedPositionOut[2]); 
    if ((cellTypes_[indexBuffer]&IS_SOLID)==0) return true; // its outside the object already 

    bool found = false; 
    const int maxIteration = 100; 

    Vector3d position; 
    Vector3d closestPoint;

    GetCellCenterPosition(ii,jj,kk,position[0],position[1],position[2]); 

    //std::cout << "inside object, position: " << position << std::endl;

    distanceField_->closestPoint(position, closestPoint);  // FIXME can I assume this point always lies on the cut? 
    Vector3d normal = distanceField_->gradient(position); 

    Vector3d reflectedPosition;
    Geometry::Reflection(position, normal, closestPoint, reflectedPosition); 

    // reflection until get out of the boundary
    int iteration = 0; 
    while(distanceField_->distance(reflectedPosition) <= distanceTolerance_ && iteration++ < maxIteration)
    {
        distanceField_->closestPoint(reflectedPosition, closestPoint); 
        normal = distanceField_->gradient(reflectedPosition); 
        Geometry::Reflection(reflectedPosition, normal, closestPoint, reflectedPosition);  // direct overwrite
    }

    if (distanceField_->distance(reflectedPosition) >= distanceTolerance_)
        found = true; 

    reflectedPositionOut << reflectedPosition.x, reflectedPosition.y, reflectedPosition.z; 

    //std::cout << "out here :) " << std::endl;
    //std::cout << "found reflection point outside object in " << iteration << " iterations, position: " << reflectedPosition << std::endl;

    return found; 

}


void UniformGridWithObject::Test_Reflection()
{
    typedef Eigen::aligned_allocator<Eigen::Vector3d> AllocatorVector3d; 

    int indexBuffer; 
    std::vector<Eigen::Vector3d, AllocatorVector3d> emittingPositions; 
    std::vector<Eigen::Vector3d, AllocatorVector3d> reflectedPositions;
    int N_success = 0; 
    int N_failure = 0; 

    std::vector<Eigen::Vector3d, AllocatorVector3d> failureReflections; 

    srand(time(NULL)); 

    Eigen::Vector3d vec3Buffer; 
    bool foundReflection; 

    // test in x direction
    for (int kk=0; kk<cellCount_[2]; kk++) 
        for (int jj=0; jj<cellCount_[1]; jj++) 
            for (int ii=0; ii<cellCount_[0]; ii++) 
            {

                FlattenIndicies(ii,jj,kk,indexBuffer); 
                if ((cellTypes_[indexBuffer] & IS_INTERFACE))
                {
                    GetCellCenterPosition(ii,jj,kk,vec3Buffer); 
                    emittingPositions.push_back(vec3Buffer); 

                    foundReflection = FlattenIndiciesWithReflection(ii,jj,kk,vec3Buffer); 
                    //GetCellCenterPosition(indexBuffer,vec3Buffer); 

                    reflectedPositions.push_back(vec3Buffer); 

                    if (foundReflection) 
                        N_success ++; 
                    else
                    {
                        N_failure ++; 
                        //failureReflections.push_back(vec3Buffer); 
                        failureReflections.push_back(emittingPositions[emittingPositions.size()-1]);
                    }

                }
            }

    std::cout << "reflection test results: \n"; 
    std::cout << " success: " << N_success << "\n"; 
    std::cout << " failure: " << N_failure << "\n"; 

    const int N_tested = N_success + N_failure; 
    Eigen::MatrixXd allEmittingPositions(N_tested,3); 
    Eigen::MatrixXd allReflectedPositions(N_tested,3); 
    Eigen::MatrixXd allColors = Eigen::MatrixXd::Random(N_tested,1); 
    Eigen::MatrixXd allFailurePositions(N_failure,3); 

    Vector3d fail; 
    for (int ii=0; ii<N_failure; ii++) 
    {
       
        fail.x=failureReflections[ii](0); 
        fail.y=failureReflections[ii](1); 
        fail.z=failureReflections[ii](2); 
        std::cout << "failure point " << ii << " has distance : " << distanceField_->distance(fail) << std::endl;
        allFailurePositions.row(ii) = failureReflections[ii]; 
    }

    for (int ii=0; ii<N_tested; ii++) 
    {
        allEmittingPositions.row(ii) = emittingPositions[ii]; 
        allReflectedPositions.row(ii) = reflectedPositions[ii]; 
    }

    allEmittingPositions.conservativeResize(Eigen::NoChange, 4); 
    allReflectedPositions.conservativeResize(Eigen::NoChange,4); 
    allEmittingPositions.col(3) = allColors; 
    allReflectedPositions.col(3) = allColors; 

    Eigen::MatrixXd allAveraged = (allEmittingPositions+allReflectedPositions)/2.0; 

    IO::writeMatrixX<double>(allFailurePositions,  "failurePositions.txt" , IO::ASCII); 
    IO::writeMatrixX<double>(allEmittingPositions,  "emittingPositions.txt" , IO::ASCII); 
    IO::writeMatrixX<double>(allReflectedPositions, "reflectedPositions.txt", IO::ASCII); 
    IO::writeMatrixX<double>(allAveraged, "averaged.txt", IO::ASCII);


}


void UniformGridWithObject::WriteCellTypes(const std::string &filename, const int &verbosity)
{

    std::shared_ptr<Eigen::MatrixXd> isSolid(new Eigen::MatrixXd(N_cells(),1)); 
    isSolid->setZero(); 

    std::shared_ptr<Eigen::MatrixXd> isInterface(new Eigen::MatrixXd(N_cells(),1)); 
    isInterface->setZero(); 

    std::shared_ptr<Eigen::MatrixXd> xInterface(new Eigen::MatrixXd(N_cells(),1)); 
    xInterface->setZero(); 

    std::shared_ptr<Eigen::MatrixXd> yInterface(new Eigen::MatrixXd(N_cells(),1)); 
    yInterface->setZero(); 

    std::shared_ptr<Eigen::MatrixXd> zInterface(new Eigen::MatrixXd(N_cells(),1)); 
    zInterface->setZero(); 

    std::shared_ptr<Eigen::MatrixXd> distanceField(new Eigen::MatrixXd(N_cells(),1)); 

    Eigen::Vector3d positionBuffer; 
    int indBuffer=0; 
    for (int kk=0; kk<cellCount_[2]; kk++) 
        for (int jj=0; jj<cellCount_[1]; jj++) 
            for (int ii=0; ii<cellCount_[0]; ii++) 
            {
                FlattenIndicies(ii,jj,kk,indBuffer); 
                GetCellCenterPosition(ii,jj,kk,positionBuffer(0),positionBuffer(1),positionBuffer(2)); 

                if (cellTypes_[indBuffer] & IS_SOLID    ) (*isSolid    )(indBuffer,0) = 1.0;
                if (cellTypes_[indBuffer] & IS_INTERFACE) (*isInterface)(indBuffer,0) = 1.0; 

                if      ((cellTypes_[indBuffer] & X_SOLID_ON_LEFT )==0 && (cellTypes_[indBuffer] & IS_INTERFACE)) (*xInterface)(indBuffer,0) =-1.0; 
                else if ((cellTypes_[indBuffer] & X_SOLID_ON_RIGHT)==0 && (cellTypes_[indBuffer] & IS_INTERFACE)) (*xInterface)(indBuffer,0) =+1.0; 

                if      ((cellTypes_[indBuffer] & Y_SOLID_ON_LEFT )==0 && (cellTypes_[indBuffer] & IS_INTERFACE)) (*yInterface)(indBuffer,0) =-1.0; 
                else if ((cellTypes_[indBuffer] & Y_SOLID_ON_RIGHT)==0 && (cellTypes_[indBuffer] & IS_INTERFACE)) (*yInterface)(indBuffer,0) =+1.0; 

                if      ((cellTypes_[indBuffer] & Z_SOLID_ON_LEFT )==0 && (cellTypes_[indBuffer] & IS_INTERFACE)) (*zInterface)(indBuffer,0) =-1.0; 
                else if ((cellTypes_[indBuffer] & Z_SOLID_ON_RIGHT)==0 && (cellTypes_[indBuffer] & IS_INTERFACE)) (*zInterface)(indBuffer,0) =+1.0; 


                (*distanceField)(indBuffer,0) = distanceField_->distance(Vector3d(positionBuffer(0),positionBuffer(1),positionBuffer(2))); 


            }

    InsertCellCenteredData("Solid", isSolid); 
    InsertCellCenteredData("Interface", isInterface); 
    InsertCellCenteredData("xInterface", xInterface); 
    InsertCellCenteredData("yInterface", yInterface); 
    InsertCellCenteredData("zInterface", zInterface); 
    InsertCellCenteredData("distanceField", distanceField); 
    if (verbosity>=0) 
    {
        WriteVTKCellCentered(filename+"_Solid", "Solid", "Solid_cells"); 
        if (verbosity>=1) 
        {
            WriteVTKCellCentered(filename+"_Interface", "Interface", "Interface_cells"); 
            WriteVTKCellCentered(filename+"_Distancefield", "distanceField", "distanceField"); 
            if (verbosity>=2)
            {
                WriteVTKCellCentered(filename+"_Interface_1", "xInterface", "Interface_components"); 
                WriteVTKCellCentered(filename+"_Interface_2", "yInterface", "Interface_components"); 
                WriteVTKCellCentered(filename+"_Interface_3", "zInterface", "Interface_components"); 
            }
        }
    }

    DeleteCellCenteredData("Solid"); 
    DeleteCellCenteredData("Interface"); 
    DeleteCellCenteredData("xInterface"); 
    DeleteCellCenteredData("yInterface"); 
    DeleteCellCenteredData("zInterface"); 
    DeleteCellCenteredData("distanceField"); 
}

void UniformGridWithObject::CellCenteredScalarHessian(const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &hessian)
{
    //std::cout << "hessian computation" << std::endl;
    //const GridData & data = GetCellCenteredData( dataName ); 
    //const int dataDimension =  data.NData(); 
    //const int NCell = data.NCells(); 
    //assert( dataDimension == 1 ); 

    //const int Nx = cellCount_[0]; 
    //const int Ny = cellCount_[1]; 
    //const int Nz = cellCount_[2]; 

    //const double dx = ( maxBound_[0] - minBound_[0] )/(double)Nx; 
    //const double dy = ( maxBound_[1] - minBound_[1] )/(double)Ny; 
    //const double dz = ( maxBound_[2] - minBound_[2] )/(double)Nz; 

    //const double dx2 = dx*dx; 
    //const double dy2 = dy*dy; 
    //const double dz2 = dz*dz; 


    //const int N_hessianNeeded = 6;
    //// determine whether to resize the input vector
    //bool resizeNeeded = false; 
    //for (auto & m : hessian) 
    //{
    //    if (nullptr == m) 
    //    {
    //        resizeNeeded = true; 
    //        break; 
    //    }

    //    if (m->rows()!=NCell || m->cols()!=1)
    //    {
    //        resizeNeeded = true; 
    //        break; 
    //    }
    //}
    //if (hessian.size()!=N_hessianNeeded)
    //    resizeNeeded = true; 

    //if (resizeNeeded) 
    //{
    //    throw std::runtime_error("**ERROR** pass in a properly resized buffer to save time");
    //    hessian.clear(); 
    //    hessian.resize(N_hessianNeeded); 
    //    for (size_t ii=0; ii<N_hessianNeeded; ii++) 
    //        hessian[ii].reset(new Eigen::MatrixXd(NCell,1));
    //}

    //int flattenedInd;
    //std::array<double,N_hessianNeeded> hessianCell; //xx, xy, xz, yy, yz, zz

    //// loop through each direction and compute the data difference and
    //// take care of the boundaries
    //#pragma omp parallel for
    //for (int kk=0; kk<Nz; kk++)
    //{
    //    for (int jj=0; jj<Ny; jj++)
    //    {
    //        for (int ii=0; ii<Ny; ii++)
    //        { 
    //            FlattenIndicies(ii,jj,kk,flattenedInd); 

    //            if (cellTypes_[flattenedInd] & IS_SOLID) continue; // do nothing

    //            hessianCell[0] = (data.Value(ii+1,jj,kk)(0) + data.Value(ii-1,jj,kk)(0) - 2.0*data.Value(ii,jj,kk)(0))/dx2; 
    //            hessianCell[3] = (data.Value(ii,jj+1,kk)(0) + data.Value(ii,jj-1,kk)(0) - 2.0*data.Value(ii,jj,kk)(0))/dy2; 
    //            hessianCell[5] = (data.Value(ii,jj,kk+1)(0) + data.Value(ii,jj,kk-1)(0) - 2.0*data.Value(ii,jj,kk)(0))/dx2; 

    //            hessianCell[1] = (data.Value(ii+1,jj+1,kk)(0) + data.Value(ii-1,jj-1,kk)(0) - data.Value(ii+1,jj-1,kk)(0) -  data.Value(ii-1,jj+1,kk)(0))/dx/dy/4.0;
    //            hessianCell[2] = (data.Value(ii+1,jj,kk+1)(0) + data.Value(ii-1,jj,kk-1)(0) - data.Value(ii+1,jj,kk-1)(0) -  data.Value(ii-1,jj,kk+1)(0))/dx/dz/4.0;
    //            hessianCell[4] = (data.Value(ii,jj+1,kk+1)(0) + data.Value(ii,jj-1,kk-1)(0) - data.Value(ii,jj-1,kk+1)(0) -  data.Value(ii,jj+1,kk-1)(0))/dy/dz/4.0;


    //            // FIXME FIXME finish it with boundary treatment.





    //                int xSign, ySign, zSign; // for boundary 

    //                if      (ii==0   )  {xSign =  1;} 
    //                else if (ii==Nx-1)  {xSign = -1;} 
    //                else                
    //                {
    //                    if      (isBoundary_[FlattenIndicies(ii-1,jj,kk)]==SOLID) xSign =  1;
    //                    else if (isBoundary_[FlattenIndicies(ii+1,jj,kk)]==SOLID) xSign = -1;
    //                    else    xSign =  0; // not on boundary 
    //                }

    //                if      (jj==0)     {ySign =  1;} 
    //                else if (jj==Ny-1)  {ySign = -1;} 
    //                else                
    //                {
    //                    if      (isBoundary_[FlattenIndicies(ii,jj-1,kk)]==SOLID) ySign =  1;
    //                    else if (isBoundary_[FlattenIndicies(ii,jj+1,kk)]==SOLID) ySign = -1;
    //                    else    ySign =  0; // not on boundary 
    //                }

    //                if      (kk==0)     {zSign =  1;} 
    //                else if (kk==Nz-1)  {zSign = -1;} 
    //                else               
    //                {
    //                    if      (isBoundary_[FlattenIndicies(ii,jj,kk-1)]==SOLID) zSign =  1;
    //                    else if (isBoundary_[FlattenIndicies(ii,jj,kk+1)]==SOLID) zSign = -1;
    //                    else    zSign =  0; // not on boundary 
    //                }

    //                // second order accurate, but it creates a little jagged
    //                // edges since it extended one more stencil to maintain 2nd
    //                // accuracy
    //                if (xSign != 0) 
    //                    buffer0 = static_cast<double>(xSign)*(-data.Value(ii,jj,kk)(0) + data.Value(ii+xSign*1,jj,kk)(0))/dx; 
    //                //buffer0 = static_cast<double>(xSign)*(-3.*data.Value(ii,jj,kk)(0) + 4.*data.Value(ii+xSign*1,jj,kk)(0) - data.Value(ii+xSign*2,jj,kk)(0))/2./dx; 
    //                else 
    //                    buffer0 = (data.Value(ii+1,jj,kk)(0) - data.Value(ii-1,jj,kk)(0))/2./dx; 

    //                if (ySign != 0) 
    //                    buffer1 = static_cast<double>(ySign)*(-data.Value(ii,jj,kk)(0) + data.Value(ii,jj+ySign*1,kk)(0))/dy; 
    //                //buffer1 = static_cast<double>(ySign)*(-3.*data.Value(ii,jj,kk)(0) + 4.*data.Value(ii,jj+ySign*1,kk)(0) - data.Value(ii,jj+ySign*2,kk)(0))/2./dy; 
    //                else 
    //                    buffer1 = (data.Value(ii,jj+1,kk)(0) - data.Value(ii,jj-1,kk)(0))/2./dy; 

    //                if (zSign != 0) 
    //                    buffer2 = static_cast<double>(zSign)*(-data.Value(ii,jj,kk)(0) + data.Value(ii,jj,kk+zSign*1)(0))/dz; 
    //                //buffer2 = static_cast<double>(zSign)*(-3.*data.Value(ii,jj,kk)(0) + 4.*data.Value(ii,jj,kk+zSign*1)(0) - data.Value(ii,jj,kk+zSign*2)(0))/2./dz; 
    //                else 
    //                    buffer2 = (data.Value(ii,jj,kk+1)(0) - data.Value(ii,jj,kk-1)(0))/2./dz; 


//#pragma omp critical
    //            {
    //                switch (i_component)
    //                {
    //                    case 0: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer0; 
    //                        (*(gradientData[1]))(flattenedInd, 0) = buffer1; 
    //                        (*(gradientData[2]))(flattenedInd, 0) = buffer2; 
    //                        break;
    //                    case 1: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer0; 
    //                        break;
    //                    case 2: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer1; 
    //                        break;
    //                    case 3: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer2; 
    //                        break;
    //                    case 12: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer0; 
    //                        (*(gradientData[1]))(flattenedInd, 0) = buffer1; 
    //                        break;
    //                    case 23: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer1; 
    //                        (*(gradientData[1]))(flattenedInd, 0) = buffer2; 
    //                        break;
    //                    case 13: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer0; 
    //                        (*(gradientData[1]))(flattenedInd, 0) = buffer2; 
    //                        break;
    //                }
    //            }
    //        } 
    //    } 
    //}

    //// local smoothing on interfacial cells to lower the impact on stencil
    //// choice on the boundaries
    //std::cout << "interface smoothing" << std::endl;
    //Eigen::VectorXd filter = SIGNAL_PROCESSING::DiscreteGaussian1D(3,1); 

    //std::vector<bool> filterMask(NCell, 1); 
    //for (int ii=0; ii<NCell; ii++) 
    //    filterMask[ii] = ((interfacialCellTypes_[ii]&IS_INTERFACE) ? true : false); 

    //for (size_t ii=0; ii<N_gradientNeeded; ii++) 
    //{
    //    InsertCellCenteredData("tmp_g", gradientData[ii]); 
    //    *(gradientData[ii]) = CellCenteredSmoothing("tmp_g", filter, filterMask); 
    //    DeleteCellCenteredData("tmp_g"); 
    //}

}


// use first order finite-difference on boundaries 
void UniformGridWithObject::CellCenteredDataGradient( const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &gradientData, const CELL_GRADIENT_COMPONENT &component)
{

    //throw std::runtime_error( "**ERROR** not implemented" ); 

    //const GridData & data = this->GetCellCenteredData( dataName ); 
    //const int dataDimension =  data.NData(); 
    //const int NCell = data.NCells(); 
    //assert( dataDimension == 1 ); 

    //const int Nx = this->cellCount_[0]; 
    //const int Ny = this->cellCount_[1]; 
    //const int Nz = this->cellCount_[2]; 

    //const double dx = ( this->maxBound_[0] - this->minBound_[0] )/(double)Nx; 
    //const double dy = ( this->maxBound_[1] - this->minBound_[1] )/(double)Ny; 
    //const double dz = ( this->maxBound_[2] - this->minBound_[2] )/(double)Nz; 

    //size_t N_gradientNeeded; 
    //const int i_component = static_cast<int>(component); 
    //if (i_component==0)
    //    N_gradientNeeded = 3; 
    //else if (i_component>0 && i_component<10)
    //    N_gradientNeeded = 1; 
    //else if (i_component>10)
    //    N_gradientNeeded = 2; 
    //else 
    //    throw std::runtime_error("**ERROR** wrong cell ponent input : "+std::to_string(i_component)); 

    //// determine whether to resize the input vector
    //bool resizeNeeded = false; 
    //for (auto & m : gradientData) 
    //{
    //    if (nullptr == m) 
    //    {
    //        resizeNeeded = true; 
    //        break; 
    //    }

    //    if (m->rows()!=NCell || m->cols()!=1)
    //    {
    //        resizeNeeded = true; 
    //        break; 
    //    }
    //}
    //if (gradientData.size()!=N_gradientNeeded)
    //    resizeNeeded = true; 

    //if (resizeNeeded) 
    //{
    //    throw std::runtime_error("**ERROR** should not resize in gradient");
    //    gradientData.clear(); 
    //    gradientData.resize(N_gradientNeeded); 
    //    for (size_t ii=0; ii<N_gradientNeeded; ii++) 
    //        gradientData[ii].reset(new Eigen::MatrixXd(NCell,1));
    //}

    //// loop through each direction and compute the data difference and
    //// take care of the boundaries
    //#pragma omp parallel for
    //for (int kk=0; kk<Nz; kk++)
    //{
    //    for (int jj=0; jj<Ny; jj++)
    //    {
    //        for (int ii=0; ii<Ny; ii++)
    //        { 
    //            const int flattenedInd = data.FlattenIndicies(ii,jj,kk); 
    //            double buffer0=0, buffer1=0, buffer2=0; 

    //            // do the computation only if its non solid
    //            if ((cellTypes_[flattenedInd] & IS_SOLID) ==0)
    //            {
    //                int xSign, ySign, zSign; // for boundary 

    //                if      (ii==0   )  {xSign =  1;} 
    //                else if (ii==Nx-1)  {xSign = -1;} 
    //                else                
    //                {
    //                    if      (cellTypes_[FlattenIndicies(ii-1,jj,kk)] & IS_SOLID) xSign =  1;
    //                    else if (cellTypes_[FlattenIndicies(ii+1,jj,kk)] & IS_SOLID) xSign = -1;
    //                    else    xSign =  0; // not on boundary 
    //                }

    //                if      (jj==0)     {ySign =  1;} 
    //                else if (jj==Ny-1)  {ySign = -1;} 
    //                else                
    //                {
    //                    if      (cellTypes_[FlattenIndicies(ii,jj-1,kk)]==SOLID) ySign =  1;
    //                    else if (cellTypes_[FlattenIndicies(ii,jj+1,kk)]==SOLID) ySign = -1;
    //                    else    ySign =  0; // not on boundary 
    //                }

    //                if      (kk==0)     {zSign =  1;} 
    //                else if (kk==Nz-1)  {zSign = -1;} 
    //                else               
    //                {
    //                    if      (cellTypes_[FlattenIndicies(ii,jj,kk-1)]==SOLID) zSign =  1;
    //                    else if (cellTypes_[FlattenIndicies(ii,jj,kk+1)]==SOLID) zSign = -1;
    //                    else    zSign =  0; // not on boundary 
    //                }

    //                // second order accurate, but it creates a little jagged
    //                // edges since it extended one more stencil to maintain 2nd
    //                // accuracy
    //                if (xSign != 0) 
    //                    buffer0 = static_cast<double>(xSign)*(-data.Value(ii,jj,kk)(0) + data.Value(ii+xSign*1,jj,kk)(0))/dx; 
    //                    //buffer0 = static_cast<double>(xSign)*(-3.*data.Value(ii,jj,kk)(0) + 4.*data.Value(ii+xSign*1,jj,kk)(0) - data.Value(ii+xSign*2,jj,kk)(0))/2./dx; 
    //                else 
    //                    buffer0 = (data.Value(ii+1,jj,kk)(0) - data.Value(ii-1,jj,kk)(0))/2./dx; 

    //                if (ySign != 0) 
    //                    buffer1 = static_cast<double>(ySign)*(-data.Value(ii,jj,kk)(0) + data.Value(ii,jj+ySign*1,kk)(0))/dy; 
    //                    //buffer1 = static_cast<double>(ySign)*(-3.*data.Value(ii,jj,kk)(0) + 4.*data.Value(ii,jj+ySign*1,kk)(0) - data.Value(ii,jj+ySign*2,kk)(0))/2./dy; 
    //                else 
    //                    buffer1 = (data.Value(ii,jj+1,kk)(0) - data.Value(ii,jj-1,kk)(0))/2./dy; 

    //                if (zSign != 0) 
    //                    buffer2 = static_cast<double>(zSign)*(-data.Value(ii,jj,kk)(0) + data.Value(ii,jj,kk+zSign*1)(0))/dz; 
    //                    //buffer2 = static_cast<double>(zSign)*(-3.*data.Value(ii,jj,kk)(0) + 4.*data.Value(ii,jj,kk+zSign*1)(0) - data.Value(ii,jj,kk+zSign*2)(0))/2./dz; 
    //                else 
    //                    buffer2 = (data.Value(ii,jj,kk+1)(0) - data.Value(ii,jj,kk-1)(0))/2./dz; 
    //            }


    //#pragma omp critical
    //            {
    //                switch (i_component)
    //                {
    //                    case 0: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer0; 
    //                        (*(gradientData[1]))(flattenedInd, 0) = buffer1; 
    //                        (*(gradientData[2]))(flattenedInd, 0) = buffer2; 
    //                        break;
    //                    case 1: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer0; 
    //                        break;
    //                    case 2: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer1; 
    //                        break;
    //                    case 3: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer2; 
    //                        break;
    //                    case 12: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer0; 
    //                        (*(gradientData[1]))(flattenedInd, 0) = buffer1; 
    //                        break;
    //                    case 23: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer1; 
    //                        (*(gradientData[1]))(flattenedInd, 0) = buffer2; 
    //                        break;
    //                    case 13: 
    //                        (*(gradientData[0]))(flattenedInd, 0) = buffer0; 
    //                        (*(gradientData[1]))(flattenedInd, 0) = buffer2; 
    //                        break;
    //                }
    //            }
    //        } 
    //    } 
    //}

    //// local smoothing on interfacial cells to lower the impact on stencil
    //// choice on the boundaries
    //std::cout << "interface smoothing" << std::endl;
    //Eigen::VectorXd filter = SIGNAL_PROCESSING::DiscreteGaussian1D(3,1); 

    //std::vector<bool> filterMask(NCell, 1); 
    //for (int ii=0; ii<NCell; ii++) 
    //    filterMask[ii] = ((interfacialCellTypes_[ii]&IS_INTERFACE) ? true : false); 

    //for (size_t ii=0; ii<N_gradientNeeded; ii++) 
    //{
    //    InsertCellCenteredData("tmp_g", gradientData[ii]); 
    //    *(gradientData[ii]) = CellCenteredSmoothing("tmp_g", filter, filterMask); 
    //    DeleteCellCenteredData("tmp_g"); 
    //}




}

