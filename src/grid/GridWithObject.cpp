#include "Grid.h"
#include "GridWithObject.h"
#include "signalprocessing/FilterDesign.h" 
#include "geometry/AffineTransformation.h" 
#include "geometry/BoundingBox.h" 
#include "distancefield/trilinearInterpolation.h"
#include "IO/IO.h"
#include "utils/STL_Wrapper.h" 
#include "math/LeastSquareSurface.h"
#include "wavesolver/PML_WaveSolver.h"
#include "wavesolver/FDTD_RigidObject.h"

template<typename T> 
void SetEigenVector3(const T &x, const T &y, const T &z, Eigen::Matrix<T,3,1> &vec3)
{
    vec3[0] = x; 
    vec3[1] = y; 
    vec3[2] = z; 
}

// recover the geometry parameters from the wave solver config file 
void UniformGridWithObject::Reinitialize(const std::string &configFile)
{
    parser_ = std::make_shared<ImpulseResponseParser>(configFile); 
    std::shared_ptr<PML_WaveSolver_Settings> settings = std::make_shared<PML_WaveSolver_Settings>(); 

    parser_->GetSolverSettings(settings); 
    parser_->GetObjects(objects_); 

    PML_WaveSolver solver(*settings, objects_); 
    BoundingBox solverBBox = solver.GetGrid().PressureBoundingBox(); 
    minBound_ = Eigen::Vector3d(solverBBox.minBound().x, solverBBox.minBound().y, solverBBox.minBound().z);
    maxBound_ = minBound_ + Eigen::Vector3d::Ones()*settings->cellSize*(double)settings->cellDivisions; 
    cellCount_[0] = settings->cellDivisions; 
    cellCount_[1] = settings->cellDivisions; 
    cellCount_[2] = settings->cellDivisions; 
    if (objects_->N()>0) 
    {
        distanceField_ = objects_->GetPtr(0)->GetSignedDistanceFieldPtr(); 
        mesh_ = objects_->GetPtr(0)->GetMeshPtr(); 
    }




    // TODO TODO 
















//    std::cout << "config file : " << configFile << std::endl;
//
//    // parse config file to get geometry parameters 
//    parser_.reset(Parser::buildParser( configFile )); 
//    if ( !parser_ ) throw std::runtime_error("**ERROR** Could not build parser from "+configFile);
//
//    mesh_.reset(parser_->getMesh("impulse_response")); 
//    if ( !mesh_ ) throw std::runtime_error("**ERROR** Could not build mesh");
//
//    solverParameters_ = parser_->getImpulseResponseParms();
//
//    double cellSize;
//    BoundingBox fieldBBox; 
//    BoundingBox solverBBox; 
//    cellCount_[0] = solverParameters_._gridResolution; 
//    cellCount_[1] = solverParameters_._gridResolution; 
//    cellCount_[2] = solverParameters_._gridResolution; 
//
//    distanceField_.reset(
//            DistanceFieldBuilder::BuildSignedClosestPointField(
//                parser_->getMeshFileName().c_str(), 
//                solverParameters_._sdfResolution, 
//                solverParameters_._sdfFilePrefix.c_str()
//                )
//            ); 
//    fieldBBox = BoundingBox(distanceField_->bmin(), distanceField_->bmax()); 
//
//    if (solverParameters_._cellSize >= 1E-14) 
//    {
//        cellSize = solverParameters_._cellSize; 
//        const REAL fieldLength = cellSize*(REAL)solverParameters_._gridResolution; 
//        const Vec3d fieldMin(-fieldLength/2.0, -fieldLength/2.0, -fieldLength/2.0); 
//        const Vec3d fieldMax(+fieldLength/2.0, +fieldLength/2.0, +fieldLength/2.0); 
//        solverBBox = BoundingBox(fieldMin, fieldMax); 
//
//    }
//    else 
//    {
//        //distanceField_.reset(
//        //        DistanceFieldBuilder::BuildSignedClosestPointField(
//        //            parser_->getMeshFileName().c_str(), 
//        //            solverParameters_._sdfResolution, 
//        //            solverParameters_._sdfFilePrefix.c_str()
//        //            )
//        //        ); 
//        //fieldBBox = BoundingBox(distanceField_->bmin(), distanceField_->bmax()); 
//        solverBBox  = BoundingBox(distanceField_->bmin(),distanceField_->bmax()); 
//        solverBBox *= solverParameters_._gridScale; 
//
//        cellSize = solverBBox.minlength()/(REAL)solverParameters_._gridResolution; 
//    }
//
//    if (!distanceField_) throw std::runtime_error("**ERROR** Could not construct distance field"); 
//
//    //fieldBBox *= solverParameters_._gridScale; 
//
//
//    minBound_ = Eigen::Vector3d(solverBBox.minBound().x, solverBBox.minBound().y, solverBBox.minBound().z);
//    maxBound_ = minBound_ + Eigen::Vector3d::Ones()*cellSize*(double)solverParameters_._gridResolution; 

    RecomputeCachedField(); 

    // reinitialize fields 

    cellTypes_.clear(); 
    cellTypes_.resize(N_cells(),0); 

    finiteDifferenceStencils_.clear(); 
    finiteDifferenceStencils_.resize(N_cells(),-1); 

    ClassifyCells(); 
    //ComputeFiniteDifferenceStencils(); 

    // output 
    std::cout << "\n" << *this << std::endl;

    //initialized_ = true; 
}

void UniformGridWithObject::ClassifyCells()
{

    std::cout << "classifying cells" << std::endl;
    const Eigen::Vector3i &cellCount = this->cellCount_; 
    const int N_cells = this->N_cells(); 


    // first rasterize the geometry
    int N_solids = 0; 
    Vector3d position; 
    for (int kk=1; kk<cellCount[2]-1; kk++) 
    {
        std::cout << " classify bulk cells. progress: " << kk << "/" << cellCount[2]-1 << "\r" << std::flush;
        for (int jj=1; jj<cellCount[1]-1; jj++) 
        {
            for (int ii=1; ii<cellCount[0]-1; ii++) 
            {
                GetCellCenterPosition(ii,jj,kk,position.x,position.y,position.z); 
                if (objects_->LowestObjectDistance(position) <= distanceTolerance_)
                {
                    cellTypes_[FlattenIndicies(ii,jj,kk)] |= IS_SOLID; 

                    // also tag the neighboring cells (this will tag a lot of
                    // cells repeatedly, but this function is executed only
                    // once
                    cellTypes_[FlattenIndicies(ii+1,jj,kk)] |= X_SOLID_ON_LEFT; 
                    cellTypes_[FlattenIndicies(ii-1,jj,kk)] |= X_SOLID_ON_RIGHT; 
                    cellTypes_[FlattenIndicies(ii,jj+1,kk)] |= Y_SOLID_ON_LEFT; 
                    cellTypes_[FlattenIndicies(ii,jj-1,kk)] |= Y_SOLID_ON_RIGHT; 
                    cellTypes_[FlattenIndicies(ii,jj,kk+1)] |= Z_SOLID_ON_LEFT; 
                    cellTypes_[FlattenIndicies(ii,jj,kk-1)] |= Z_SOLID_ON_RIGHT; 
                    N_solids ++; 
                }
            }
        }
    }
    std::cout << std::endl;


    // then find the interfacial cells: solid cell who has neighbors of mix
    // types    
    int index; 
    int N_interfacial = 0; 
    int countBuffer;
    for (int kk=1; kk<cellCount[2]-1; kk++) 
    {
        std::cout << " classify interfacial cells. progress: " << kk << "/" << cellCount[2]-2 << "\r" << std::flush;
        for (int jj=1; jj<cellCount[1]-1; jj++) 
        {
            for (int ii=1; ii<cellCount[0]-1; ii++) 
            {
                index = FlattenIndicies(ii,jj,kk); 

                countBuffer = 0; 
                if (cellTypes_[index] & X_SOLID_ON_LEFT)  countBuffer+=1; 
                else                                      countBuffer-=1; 
                if (cellTypes_[index] & X_SOLID_ON_RIGHT) countBuffer+=1; 
                else                                      countBuffer-=1; 
                if (cellTypes_[index] & Y_SOLID_ON_LEFT)  countBuffer+=1; 
                else                                      countBuffer-=1; 
                if (cellTypes_[index] & Y_SOLID_ON_RIGHT) countBuffer+=1; 
                else                                      countBuffer-=1; 
                if (cellTypes_[index] & Z_SOLID_ON_LEFT)  countBuffer+=1; 
                else                                      countBuffer-=1; 
                if (cellTypes_[index] & Z_SOLID_ON_RIGHT) countBuffer+=1; 
                else                                      countBuffer-=1; 

                if (countBuffer!=6 && countBuffer!=-6) 
                {
                    cellTypes_[index] |= IS_INTERFACE; 
                    N_interfacial++; 
                }


                //if (cellTypes_[FlattenIndicies(ii-1,jj  ,kk  )] & IS_SOLID) cellTypes_[index] |= X_SOLID_ON_LEFT ; 
                //if (cellTypes_[FlattenIndicies(ii+1,jj  ,kk  )] & IS_SOLID) cellTypes_[index] |= X_SOLID_ON_RIGHT;
                //if (cellTypes_[FlattenIndicies(ii  ,jj-1,kk  )] & IS_SOLID) cellTypes_[index] |= Y_SOLID_ON_LEFT ; 
                //if (cellTypes_[FlattenIndicies(ii  ,jj+1,kk  )] & IS_SOLID) cellTypes_[index] |= Y_SOLID_ON_RIGHT; 
                //if (cellTypes_[FlattenIndicies(ii  ,jj  ,kk-1)] & IS_SOLID) cellTypes_[index] |= Z_SOLID_ON_LEFT ; 
                //if (cellTypes_[FlattenIndicies(ii  ,jj  ,kk+1)] & IS_SOLID) cellTypes_[index] |= Z_SOLID_ON_RIGHT; 
                //if ((cellTypes_[index] & X_SOLID_ON_LEFT)==0 || (cellTypes_[index] & X_SOLID_ON_RIGHT)==0 || 
                //    (cellTypes_[index] & Y_SOLID_ON_LEFT)==0 || (cellTypes_[index] & Y_SOLID_ON_RIGHT)==0 || 
                //    (cellTypes_[index] & Z_SOLID_ON_LEFT)==0 || (cellTypes_[index] & Z_SOLID_ON_RIGHT)==0)
                //{
                //    cellTypes_[index] |= IS_INTERFACE; 
                //    N_interfacial++; 
                //}

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
bool UniformGridWithObject::FlattenIndiciesWithReflection(const int &ii, const int &jj, const int &kk, int &nearestCellIndex, Eigen::Vector3d &reflectedPositionOut) const
{

    int indexBuffer; 
    FlattenIndicies(ii,jj,kk,indexBuffer); 
    nearestCellIndex = indexBuffer; // return the original cell index if anything is wrong
    GetCellCenterPosition(ii,jj,kk,reflectedPositionOut[0],reflectedPositionOut[1],reflectedPositionOut[2]); 
    if ((cellTypes_[indexBuffer]&IS_SOLID)==0) return true; // its outside the object already 

    const int maxIteration = 2;

    Vector3d position; 
    Vector3d closestPoint;
    Vector3d normal;

    GetCellCenterPosition(ii,jj,kk,position[0],position[1],position[2]); 

    //std::cout << "inside object, position: " << position << std::endl;

    assert(objects_->N()>0); 
    FDTD_RigidObject &object = objects_->Get(0); 
    double distance = object.DistanceToMesh(position); 
    object.NormalToMesh(position, normal); 
    normal.normalize(); 
    closestPoint = position - normal*distance; 
    //distanceField_->closestPoint(position, closestPoint);  // FIXME can I assume this point always lies on the cut? 
    //Vector3d normal = distanceField_->gradient(position); 

    Vector3d reflectedPosition;
    Geometry::Reflection(position, normal, closestPoint, reflectedPosition); 

    // reflection until get out of the boundary
    int iteration = 0; 
    while(object.DistanceToMesh(reflectedPosition) <= distanceTolerance_ && iteration++ < maxIteration)
    {
        distance = object.DistanceToMesh(reflectedPosition); 
        object.NormalToMesh(reflectedPosition, normal); 
        normal.normalize(); 
        closestPoint = position - normal*distance; 
        Geometry::Reflection(reflectedPosition, normal, closestPoint, reflectedPosition);  // direct overwrite
    }



    //std::cout << "out here :) " << std::endl;
    //std::cout << "found reflection point outside object in " << iteration << " iterations, position: " << reflectedPosition << std::endl;

   
    // cannot find it then fallback to the cell index.
    if (object.DistanceToMesh(reflectedPosition)<=distanceTolerance_)
    {
        nearestCellIndex = indexBuffer; 
        return false; 
    }
    else 
    {
        reflectedPositionOut << reflectedPosition.x, reflectedPosition.y, reflectedPosition.z; 
    }

    // find enclosing cell
    assert(reflectedPosition.x>=minBound_[0] && reflectedPosition.y>=minBound_[1] && reflectedPosition.z>=minBound_[2]); 
    const int xIndex = max(min((int)((reflectedPosition.x-minBound_[0])/dx_[0]), cellCount_[0]-1), 1);  // slightly more restrictive clamp to make BC
    const int yIndex = max(min((int)((reflectedPosition.y-minBound_[1])/dx_[1]), cellCount_[1]-1), 1); 
    const int zIndex = max(min((int)((reflectedPosition.z-minBound_[2])/dx_[2]), cellCount_[2]-1), 1); 

    // search the neighbor 5x5 to find the closest non-solid cells
    double closestDistance=std::numeric_limits<double>::max(); 
    double distanceBuffer; 
    nearestCellIndex = -1; 
    int searchProximity = 2; 
    for (int i=-searchProximity; i<searchProximity; i++) 
        for (int j=-searchProximity; j<searchProximity; j++) 
            for (int k=-searchProximity; k<searchProximity; k++) 
            {
                FlattenIndicies(xIndex+i, yIndex+j, zIndex+k, indexBuffer); 

                if ((cellTypes_[indexBuffer] & IS_SOLID)==0) 
                {
                    GetCellCenterPosition(indexBuffer, position.x, position.y, position.z); 
                    distanceBuffer = (position-reflectedPosition).length(); 
                    if (distanceBuffer<closestDistance)
                    {
                        closestDistance = distanceBuffer; 
                        nearestCellIndex = indexBuffer; 
                    }
                    
                }
                    
            }

    if (nearestCellIndex>=0)
        return true; 
    else  // fallback to cell index if in the close proximity of reflection point there is no non-solid cell
    {
        FlattenIndicies(ii,jj,kk,nearestCellIndex); 
        return false; 
    }

    return false; 
}

int UniformGridWithObject::FindKNearestFluidCells(const int &K, const Eigen::Vector3d &centerPosition, std::vector<int> &nearestNeighbors) const 
{

    const int searchProximity = 4;  // only search k cells in proximity

    const bool boundsCheckPass= CheckLowerBounds(centerPosition[0]-searchProximity*dx_[0],centerPosition[1]-searchProximity*dx_[1],centerPosition[2]-searchProximity*dx_[2]) && 
                                CheckUpperBounds(centerPosition[0]+searchProximity*dx_[0],centerPosition[1]+searchProximity*dx_[1],centerPosition[2]+searchProximity*dx_[2]); 
    if (!boundsCheckPass) 
        throw std::runtime_error("**ERROR** center position for k nearest point search out of bounds. Could be that object is too close to grid boundary. Input center position: " + to_string(centerPosition[0]) + ", " + std::to_string(centerPosition[1]) + ", " + std::to_string(centerPosition[2])); 

    // the previous check ensures this is greater than zero
    const int xIndex = (int)((centerPosition[0]-minBound_[0])/dx_[0]);
    const int yIndex = (int)((centerPosition[1]-minBound_[1])/dx_[1]); 
    const int zIndex = (int)((centerPosition[2]-minBound_[2])/dx_[2]); 

    // search the K-by-K neighbor to find the closest cells
    int indexBuffer; 
    nearestNeighbors.clear(); 
    std::vector<double> distances; 
    Eigen::Vector3d positionBuffer; 
    double distanceBuffer; 

    // TODO if it works, can use KD tree to accelerate the search 
    for (int i=-searchProximity; i<searchProximity; i++) 
        for (int j=-searchProximity; j<searchProximity; j++) 
            for (int k=-searchProximity; k<searchProximity; k++) 
            {
                FlattenIndicies(xIndex+i, yIndex+j, zIndex+k, indexBuffer); 
                if ((cellTypes_[indexBuffer] & IS_SOLID)==0)  // record both index and distance 
                {
                    GetCellCenterPosition(indexBuffer, positionBuffer[0], positionBuffer[1], positionBuffer[2]); 
                    distanceBuffer = (positionBuffer-centerPosition).squaredNorm();  // squared should be faster

                    nearestNeighbors.push_back(indexBuffer); 
                    distances.push_back(distanceBuffer); 
                }
                    
            }

    std::vector<size_t> permutations = STL_Wrapper::SortIndicies(distances); 
    std::vector<int> indiciesBuffer = nearestNeighbors; 
    for (size_t ii=0; ii<indiciesBuffer.size(); ii++)
        nearestNeighbors[ii] = indiciesBuffer[permutations[ii]]; 

    const size_t numberNearestIndiciesFound = std::min<size_t>(nearestNeighbors.size(), (size_t)K);

    if (numberNearestIndiciesFound < nearestNeighbors.size())
        nearestNeighbors.erase(nearestNeighbors.begin()+numberNearestIndiciesFound, nearestNeighbors.end()); 


    return (int)numberNearestIndiciesFound; 

}


void UniformGridWithObject::ComputeFiniteDifferenceStencils()
{

    std::cout << "computing finite-difference stencils for solid cells" << std::endl;
    int N_failure = 0; 
    int N_success = 0; 
    int count = 0; 

#pragma omp parallel for 
    for (int kk=0; kk<cellCount_[2]; kk++) 
    {
#pragma omp critical
        {
            std::cout << " progress : " << count << "/" << cellCount_[2]-1 << "\r" << std::flush;
            count ++; 
        }

        int cellIndex;
        int stencilIndex; 
        Eigen::Vector3d posBuffer;
        for (int jj=0; jj<cellCount_[1]; jj++) 
        {
            for (int ii=0; ii<cellCount_[0]; ii++) 
            {
                FlattenIndicies(ii,jj,kk,cellIndex); 
                if ((cellTypes_[cellIndex] & IS_SOLID))
                {
                    bool foundReflection = FlattenIndiciesWithReflection(ii,jj,kk,stencilIndex,posBuffer); 

                    if (foundReflection) 
                    {
                        N_success ++; 
                    }
                    else 
                    {
                        stencilIndex = cellIndex; // fall back to the cell index if reflection attempt failed 
                        N_failure ++; 
                    }
                }
                else // if its not solid, use cell index
                {
                    stencilIndex = cellIndex; 
                }

#pragma omp critical 
                finiteDifferenceStencils_[cellIndex] = stencilIndex; 
            }
        }
    }
    std::cout << std::endl;

    std::cout << " stencils cached\n"; 
    std::cout << "  success reflection : " << N_success << "\n";
    std::cout << "  failure reflection : " << N_failure << "\n";
    std::cout << "    * in the failure case, program falls back to use original cell index as stencil" << std::endl; // increase sdf octree depth or have less concave geometry would help


}


void UniformGridWithObject::Test_Reflection()
{

    std::cout << "testing reflection stencils " << std::endl;

    typedef Eigen::aligned_allocator<Eigen::Vector3d> AllocatorVector3d; 

    int indexBuffer; 
    std::vector<Eigen::Vector3d, AllocatorVector3d> emittingPositions; 
    std::vector<Eigen::Vector3d, AllocatorVector3d> reflectedPositions;
    std::vector<Eigen::Vector3d, AllocatorVector3d> nearestNeighborPositions;
    int N_success = 0; 
    int N_failure = 0; 

    std::vector<Eigen::Vector3d, AllocatorVector3d> failureReflections; 

    srand(time(NULL)); 

    Eigen::Vector3d vec3Buffer; 
    bool foundReflection; 

    const int computeLayers = 1; 

    // test in x direction
    for (int kk=computeLayers; kk<cellCount_[2]-computeLayers; kk++) 
    {
        std::cout << kk << "\r" << std::flush;
        for (int jj=computeLayers; jj<cellCount_[1]-computeLayers; jj++) 
            for (int ii=computeLayers; ii<cellCount_[0]-computeLayers; ii++) 
            {

                FlattenIndicies(ii,jj,kk,indexBuffer); 


                if (cellTypes_[indexBuffer] & IS_INTERFACE)
                {
                    int mm = 0; 
                    //for (int mm=-computeLayers; mm<computeLayers; mm++)
                    {
                        // x extension
                        //
                        FlattenIndicies(ii+mm,jj,kk,indexBuffer); 
                        if (cellTypes_[indexBuffer] & IS_SOLID)
                        {
                            GetCellCenterPosition(ii+mm,jj,kk,vec3Buffer); 
                            emittingPositions.push_back(vec3Buffer); 

                            foundReflection = FlattenIndiciesWithReflection(ii+mm,jj,kk,indexBuffer,vec3Buffer); 
                            reflectedPositions.push_back(vec3Buffer); 

                            GetCellCenterPosition(indexBuffer,vec3Buffer); 
                            nearestNeighborPositions.push_back(vec3Buffer);

                            if (foundReflection) 
                                N_success ++; 
                            else
                            {
                                N_failure ++; 
                                //failureReflections.push_back(vec3Buffer); 
                                failureReflections.push_back(emittingPositions[emittingPositions.size()-1]);
                            }
                        }

                        // y extension
                        FlattenIndicies(ii,jj+mm,kk,indexBuffer); 
                        if (cellTypes_[indexBuffer] & IS_SOLID)
                        {
                            GetCellCenterPosition(ii,jj+mm,kk,vec3Buffer); 
                            emittingPositions.push_back(vec3Buffer); 

                            foundReflection = FlattenIndiciesWithReflection(ii,jj+mm,kk,indexBuffer,vec3Buffer); 
                            reflectedPositions.push_back(vec3Buffer); 

                            GetCellCenterPosition(indexBuffer,vec3Buffer); 
                            nearestNeighborPositions.push_back(vec3Buffer); 


                            if (foundReflection) 
                                N_success ++; 
                            else
                            {
                                N_failure ++; 
                                //failureReflections.push_back(vec3Buffer); 
                                failureReflections.push_back(emittingPositions[emittingPositions.size()-1]);
                            }
                        }

                        // z extension
                        FlattenIndicies(ii,jj,kk+mm,indexBuffer); 
                        if (cellTypes_[indexBuffer] & IS_SOLID)
                        {
                            GetCellCenterPosition(ii,jj,kk+mm,vec3Buffer); 
                            emittingPositions.push_back(vec3Buffer); 

                            foundReflection = FlattenIndiciesWithReflection(ii,jj,kk+mm,indexBuffer,vec3Buffer); 
                            reflectedPositions.push_back(vec3Buffer); 

                            GetCellCenterPosition(indexBuffer,vec3Buffer); 
                            nearestNeighborPositions.push_back(vec3Buffer); 


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
                }

            }
    }
    std::cout << std::endl;

    std::cout << "reflection test results: \n"; 
    std::cout << " success: " << N_success << "\n"; 
    std::cout << " failure: " << N_failure << "\n"; 

    const int N_tested = N_success + N_failure; 
    Eigen::MatrixXd allEmittingPositions(N_tested,3); 
    Eigen::MatrixXd allReflectedPositions(N_tested,3); 
    Eigen::MatrixXd allNearestNeighborPositions(N_tested,3); 
    Eigen::MatrixXd allColors = Eigen::MatrixXd::Random(N_tested,1); 
    Eigen::MatrixXd allFailurePositions(N_failure,3); 

    Vector3d fail; 
    assert(objects_->N()>0); 
    FDTD_RigidObject &object = objects_->Get(0); 
    for (int ii=0; ii<N_failure; ii++) 
    {

        fail.x=failureReflections[ii](0); 
        fail.y=failureReflections[ii](1); 
        fail.z=failureReflections[ii](2); 
        std::cout << "failure point " << ii << " has distance : " << object.DistanceToMesh(fail) << std::endl;
        allFailurePositions.row(ii) = failureReflections[ii]; 
    }

    for (int ii=0; ii<N_tested; ii++) 
    {
        allEmittingPositions.row(ii) = emittingPositions[ii]; 
        allReflectedPositions.row(ii) = reflectedPositions[ii]; 
        allNearestNeighborPositions.row(ii) = nearestNeighborPositions[ii]; 
    }

    allEmittingPositions.conservativeResize(Eigen::NoChange, 4); 
    allReflectedPositions.conservativeResize(Eigen::NoChange,4); 
    allEmittingPositions.col(3) = allColors; 
    allReflectedPositions.col(3) = allColors; 

    Eigen::MatrixXd allAveraged = (allEmittingPositions+allReflectedPositions)/2.0; 

    IO::writeMatrixX<double>(allFailurePositions,  "failurePositions.txt" , IO::ASCII); 
    IO::writeMatrixX<double>(allEmittingPositions,  "emittingPositions.txt" , IO::ASCII); 
    IO::writeMatrixX<double>(allReflectedPositions, "reflectedPositions.txt", IO::ASCII); 
    IO::writeMatrixX<double>(allNearestNeighborPositions, "nearestNeighborPositions.txt", IO::ASCII); 
    IO::writeMatrixX<double>(allAveraged, "averaged.txt", IO::ASCII);


}


void UniformGridWithObject::WriteCellTypes(const std::string &filename, const int &verbosity)
{

    std::cout << "writing cell types " << std::endl;

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

    assert(objects_->N()>0); 
    FDTD_RigidObject &object = objects_->Get(0); 

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


                (*distanceField)(indBuffer,0) = object.DistanceToMesh(Vector3d(positionBuffer(0),positionBuffer(1),positionBuffer(2))); 
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

int UniformGridWithObject::GetStencilIndex(Eigen::Vector3i &indicies)
{
    if (!initialized_)
        ComputeFiniteDifferenceStencils(); 

    int buffer; 
    GetStencilIndex(indicies,buffer); 

    return buffer;
}

void UniformGridWithObject::GetStencilIndex(Eigen::Vector3i &indicies, int &stencilIndex)
{
    // clamp the stencils 
    indicies[0] = min(max(indicies[0],0),cellCount_[0]-1); 
    indicies[1] = min(max(indicies[1],0),cellCount_[1]-1); 
    indicies[2] = min(max(indicies[2],0),cellCount_[2]-1); 

    int flattenIndex = FlattenIndicies(indicies[0], indicies[1], indicies[2]); 

    stencilIndex = finiteDifferenceStencils_[flattenIndex]; 
}


double UniformGridWithObject::GetStencilDataScalar(const GridData &data, const int &ii, const int &jj, const int &kk, const ReflectionFetchMethod &reflectionFetchMethod) const 
{
    assert(data.NData()==1);


    // clamp the stencils
    const int xIndex = min(max(ii,0),cellCount_[0]-1); 
    const int yIndex = min(max(jj,0),cellCount_[1]-1); 
    const int zIndex = min(max(kk,0),cellCount_[2]-1); 

    int indexBuffer;
    FlattenIndicies(xIndex,yIndex,zIndex,indexBuffer); 

    // early return if its interior cell
    if ((cellTypes_[indexBuffer]&IS_SOLID)==0)
        return data.Value(indexBuffer)(0); 

    Eigen::Vector3d reflectedPosition; 

    // fetch the reflected cell position
    FlattenIndiciesWithReflection(xIndex,yIndex,zIndex,indexBuffer,reflectedPosition); 

    double result = 0; 

    switch (reflectionFetchMethod) 
    {
        case TRILINEAR:
        {
            // find the largest cell index that is smaller than the reflected position
            Eigen::Vector3i lowestCellIndex = ((reflectedPosition-dx_/2.0).cwiseQuotient(dx_)).cast<int>(); 
            lowestCellIndex[0] = min(max(lowestCellIndex[0],0),cellCount_[0]-2); 
            lowestCellIndex[1] = min(max(lowestCellIndex[1],0),cellCount_[1]-2); 
            lowestCellIndex[2] = min(max(lowestCellIndex[2],0),cellCount_[2]-2); 

            FlattenIndicies(lowestCellIndex[0], lowestCellIndex[1], lowestCellIndex[2], indexBuffer); 

            Eigen::Vector3d lowestCellPosition  = GetCellCenterPosition(lowestCellIndex); 
            Eigen::Vector3d highestCellPosition = GetCellCenterPosition(lowestCellIndex+Eigen::Vector3i::Ones()); 

            Eigen::Vector3d weights = (reflectedPosition-lowestCellPosition).cwiseQuotient(highestCellPosition-lowestCellPosition); 

            const Eigen::Vector3i &l = lowestCellIndex; 
            const double v000 = data.Value(l[0]  ,l[1]  ,l[2]  )(0); 
            const double v100 = data.Value(l[0]+1,l[1]  ,l[2]  )(0); 
            const double v110 = data.Value(l[0]+1,l[1]+1,l[2]  )(0); 
            const double v010 = data.Value(l[0]  ,l[1]+1,l[2]  )(0); 
            const double v001 = data.Value(l[0]  ,l[1]  ,l[2]+1)(0); 
            const double v101 = data.Value(l[0]+1,l[1]  ,l[2]+1)(0); 
            const double v111 = data.Value(l[0]+1,l[1]+1,l[2]+1)(0); 
            const double v011 = data.Value(l[0]  ,l[1]+1,l[2]+1)(0); 

            result = TRILINEAR_INTERPOLATION(weights[0],weights[1],weights[2],v000,v100,v110,v010,v001,v101,v111,v011); 
            break; 
        }


        // find four nearest neighbor closest to the reflected position and fit
        // a linear surface then evaulate the function at the reflected
        // position.
        case LEAST_SQUARE: 
        {
            const int N_samples = 12; 
            std::vector<int> foundCellIndicies; 

            const int numberFoundCells = FindKNearestFluidCells(N_samples, reflectedPosition, foundCellIndicies); 

            if (numberFoundCells < N_samples)
                std::cout << "**WARNING** number Found cells = " << numberFoundCells << std::endl;

            if (numberFoundCells < 4) throw std::runtime_error("**ERROR** number of neighbors found was less than enough to fit a linear polynomial"); 

            Eigen::MatrixXd samplePoints(numberFoundCells,3); 
            Eigen::VectorXd sampleValues(numberFoundCells); 
            for (int ii=0; ii<numberFoundCells; ii++) 
            {
                GetCellCenterPosition(foundCellIndicies[ii], samplePoints(ii,0), samplePoints(ii,1), samplePoints(ii,2)); 
                sampleValues(ii) = data.Value(foundCellIndicies[ii])(0); 
            }
            LeastSquareSurfaceLinear3D fittedSurface; 
            fittedSurface.ComputeCoefficients(samplePoints, sampleValues); 
            result = fittedSurface.Evaluate(reflectedPosition); 


            std::cout << "\n"; 
            STL_Wrapper::PrintVectorContent(std::cout, foundCellIndicies); 
            std::cout << "solid:";
            for (int ii=0; ii<numberFoundCells; ii++) 
                std::cout << (cellTypes_[foundCellIndicies[ii]]&IS_SOLID) << " "; 
            std::cout << "\n";

            //std::cout << "sample points" << samplePoints << std::endl;
            std::cout << "sample values" << sampleValues.transpose() << std::endl;
            std::cout << "result = " << result << std::endl;

            break;
        }

    }

    return result; 

}


void UniformGridWithObject::CellCenteredScalarHessian(const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &hessian)
{
    std::cout << "hessian computation" << std::endl;
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

    const double dx2 = dx*dx; 
    const double dy2 = dy*dy; 
    const double dz2 = dz*dz; 


    const int N_hessianNeeded = 6;
    // determine whether to resize the input vector
    bool resizeNeeded = false; 
    for (auto & m : hessian) 
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
    if (hessian.size()!=N_hessianNeeded)
        resizeNeeded = true; 

    if (resizeNeeded) 
    {
        throw std::runtime_error("**ERROR** please pass in a properly resized buffer to save time");
        hessian.clear(); 
        hessian.resize(N_hessianNeeded); 
        for (size_t ii=0; ii<N_hessianNeeded; ii++) 
            hessian[ii].reset(new Eigen::MatrixXd(NCell,1));
    }

    int flattenedInd;
    std::array<double,N_hessianNeeded> hessianCell; //xx, xy, xz, yy, yz, zz
    Eigen::Vector3i indiciesBuffer; 

    const ReflectionFetchMethod fetchMethod = LEAST_SQUARE; 
    //const ReflectionFetchMethod fetchMethod = TRILINEAR; 

    // loop through each direction and compute the data difference and
    // take care of the boundaries
    for (int kk=0; kk<Nz; kk++)
    {
        //std::cout << kk << "/" << Nz-1 << "\r" << std::flush;
        for (int jj=0; jj<Ny; jj++)
        {
            for (int ii=0; ii<Nx; ii++)
            { 
                FlattenIndicies(ii,jj,kk,flattenedInd); 

                if (cellTypes_[flattenedInd] & IS_SOLID) continue; // do nothing if it is solid

                //if (ii==66 && jj==50 && kk==49)
                //    std::cout << "cell position = " << GetCellCenterPosition(ii,jj,kk) << std::endl;



                //debug 
                //hessianCell[0] = GetStencilDataScalar(data,ii+1,jj,kk,TRILINEAR) - GetStencilDataScalar(data,ii-1,jj,kk,TRILINEAR); 
                //hessianCell[0] = GetStencilDataScalar(data,ii+1,jj,kk,LEAST_SQUARE) - GetStencilDataScalar(data,ii-1,jj,kk,LEAST_SQUARE); 
                //hessianCell[0] /= (2.0*dx_[0]); 

                //hessianCell[0] = 0; 
                //SetEigenVector3<int>(ii+1,jj  ,kk  , indiciesBuffer); 
                //hessianCell[0] += data.Value(GetStencilIndex(indiciesBuffer))(0); 
                //SetEigenVector3<int>(ii-1,jj  ,kk  , indiciesBuffer); 
                //hessianCell[0] -= data.Value(GetStencilIndex(indiciesBuffer))(0); 
                //hessianCell[0] /= (2.0*dx_[0]); 




                // xx
                hessianCell[0] = GetStencilDataScalar(data,ii+1,jj,kk,fetchMethod) + GetStencilDataScalar(data,ii-1,jj,kk,fetchMethod) - 2.0*GetStencilDataScalar(data,ii,jj,kk,fetchMethod); 
                hessianCell[0] /= dx2; 

                // yy
                hessianCell[3] = GetStencilDataScalar(data,ii,jj+1,kk,fetchMethod) + GetStencilDataScalar(data,ii,jj-1,kk,fetchMethod) - 2.0*GetStencilDataScalar(data,ii,jj,kk,fetchMethod); 
                hessianCell[3] /= dy2; 


                // zz
                hessianCell[5] = GetStencilDataScalar(data,ii,jj,kk+1,fetchMethod) + GetStencilDataScalar(data,ii,jj,kk-1,fetchMethod) - 2.0*GetStencilDataScalar(data,ii,jj,kk,fetchMethod); 
                hessianCell[5] /= dz2; 


                // xy
                hessianCell[1]  = GetStencilDataScalar(data,ii+1,jj+1,kk,fetchMethod) + GetStencilDataScalar(data,ii-1,jj-1,kk,fetchMethod)
                                - GetStencilDataScalar(data,ii+1,jj-1,kk,fetchMethod) - GetStencilDataScalar(data,ii-1,jj+1,kk,fetchMethod); 
                hessianCell[1] /= (dx*dy*4.0);


                // xz
                hessianCell[2]  = GetStencilDataScalar(data,ii+1,jj,kk+1,fetchMethod) + GetStencilDataScalar(data,ii-1,jj,kk-1,fetchMethod)
                                - GetStencilDataScalar(data,ii+1,jj,kk-1,fetchMethod) - GetStencilDataScalar(data,ii-1,jj,kk+1,fetchMethod); 
                hessianCell[2] /= (dx*dz*4.0);

                // yz
                hessianCell[4]  = GetStencilDataScalar(data,ii,jj+1,kk+1,fetchMethod) + GetStencilDataScalar(data,ii,jj-1,kk-1,fetchMethod)
                                - GetStencilDataScalar(data,ii,jj+1,kk-1,fetchMethod) - GetStencilDataScalar(data,ii,jj-1,kk+1,fetchMethod); 
                hessianCell[4] /= (dy*dz*4.0);

                /* using nearest neighbor
                // xx
                hessianCell[0] = 0; 
                SetEigenVector3<int>(ii+1,jj  ,kk  , indiciesBuffer); 
                hessianCell[0] += data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii-1,jj  ,kk  , indiciesBuffer); 
                hessianCell[0] += data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii  ,jj  ,kk  , indiciesBuffer); 
                hessianCell[0] += -2.0*data.Value(GetStencilIndex(indiciesBuffer))(0); 
                hessianCell[0] /= dx2; 

                // yy
                hessianCell[3] = 0; 
                SetEigenVector3<int>(ii  ,jj+1,kk  , indiciesBuffer); 
                hessianCell[3] += data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii  ,jj-1,kk  , indiciesBuffer); 
                hessianCell[3] += data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii  ,jj  ,kk  , indiciesBuffer); 
                hessianCell[3] += -2.0*data.Value(GetStencilIndex(indiciesBuffer))(0); 
                hessianCell[3] /= dy2; 

                // zz
                hessianCell[5] = 0; 
                SetEigenVector3<int>(ii  ,jj  ,kk+1, indiciesBuffer); 
                hessianCell[5] += data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii  ,jj  ,kk-1, indiciesBuffer); 
                hessianCell[5] += data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii  ,jj  ,kk  , indiciesBuffer); 
                hessianCell[5] += -2.0*data.Value(GetStencilIndex(indiciesBuffer))(0); 
                hessianCell[5] /= dz2; 

                // xy
                hessianCell[1] = 0; 
                SetEigenVector3<int>(ii+1,jj+1,kk  , indiciesBuffer); 
                hessianCell[1] +=  data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii-1,jj-1,kk  , indiciesBuffer); 
                hessianCell[1] +=  data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii+1,jj-1,kk  , indiciesBuffer); 
                hessianCell[1] += -data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii-1,jj+1,kk  , indiciesBuffer); 
                hessianCell[1] += -data.Value(GetStencilIndex(indiciesBuffer))(0); 
                hessianCell[1] /= (dx*dy*4.0); 

                // xz
                hessianCell[2] = 0; 
                SetEigenVector3<int>(ii+1,jj  ,kk+1, indiciesBuffer); 
                hessianCell[2] +=  data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii-1,jj  ,kk-1, indiciesBuffer); 
                hessianCell[2] +=  data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii+1,jj  ,kk-1, indiciesBuffer); 
                hessianCell[2] += -data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii-1,jj  ,kk+1, indiciesBuffer); 
                hessianCell[2] += -data.Value(GetStencilIndex(indiciesBuffer))(0); 
                hessianCell[2] /= (dx*dz*4.0); 

                // yz
                hessianCell[4] = 0; 
                SetEigenVector3<int>(ii  ,jj+1,kk+1, indiciesBuffer); 
                hessianCell[4] +=  data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii  ,jj-1,kk-1, indiciesBuffer); 
                hessianCell[4] +=  data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii  ,jj+1,kk-1, indiciesBuffer); 
                hessianCell[4] += -data.Value(GetStencilIndex(indiciesBuffer))(0); 
                SetEigenVector3<int>(ii  ,jj-1,kk+1, indiciesBuffer); 
                hessianCell[4] += -data.Value(GetStencilIndex(indiciesBuffer))(0); 
                hessianCell[4] /= (dy*dz*4.0); 
                */


                {
                    (*(hessian[0]))(flattenedInd, 0) = hessianCell[0]; 
                    (*(hessian[1]))(flattenedInd, 0) = hessianCell[1]; 
                    (*(hessian[2]))(flattenedInd, 0) = hessianCell[2]; 
                    (*(hessian[3]))(flattenedInd, 0) = hessianCell[3]; 
                    (*(hessian[4]))(flattenedInd, 0) = hessianCell[4]; 
                    (*(hessian[5]))(flattenedInd, 0) = hessianCell[5]; 
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
        filterMask[ii] = ((cellTypes_[ii]&IS_INTERFACE) ? true : false); 

    for (size_t ii=0; ii<hessian.size(); ii++) 
    {
        InsertCellCenteredData("tmp_h", hessian[ii]); 
        *(hessian[ii]) = CellCenteredSmoothing("tmp_h", filter, filterMask); 
        DeleteCellCenteredData("tmp_h"); 
    }
    

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

