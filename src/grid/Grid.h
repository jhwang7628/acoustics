#ifndef GRID_H
#define GRID_H 

#include <iostream>
#include <string>
#include <unordered_map>
#include <iomanip> 
#include "vtkConverter/vtkConverter.h"

#include <Eigen/Dense> 

template <class T> class GridData; 
template <class T> class Grid;
template <class T> 
Eigen::VectorXd TrilinearInterpolate( const Eigen::Vector3d& xmin, const Eigen::Vector3d& xmax, const Eigen::Vector3i& N, const Eigen::Vector3d& position, const GridData<T>& data );

/*
 * data defined on the grid. each row in the dat matrix represents data
 * associated to one cell(vertex). therefore the number of rows should be equal to 
 * the total number of cells(vertex). 
 *
 * the storage order is assume to be first x then y then z. e.g. the R3
 * indicies corresponding to a 8 rows data matrix is something like
 *
 * x y z
 * -----
 * 0 0 0
 * 1 0 0
 * 0 1 0
 * 1 1 0
 * 0 0 1
 * 1 0 1
 * 0 1 1
 * 1 1 1
 *
 * see the getter methods GetValue() for more details
 *
 *
 */
template <class T>
class GridData
{
    private: 
        std::string name;
        std::shared_ptr<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> dat; 

        /// the grid that owns this data
        Grid<T>& owner_; 

    public: 

        GridData( const std::string& nameIn, std::shared_ptr<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> dataIn, Grid<T>& owner ) : name(nameIn), dat(dataIn), owner_(owner) { } 

        /// getters 
        Eigen::Matrix<T, 1, Eigen::Dynamic> Value( const int & index ) const { return dat->row(index); }
        Eigen::Matrix<T, 1, Eigen::Dynamic> Value( const int xIndex, const int yIndex, const int zIndex ) const 
        { 
            const int ind = FlattenIndicies(xIndex, yIndex, zIndex);
            return dat->row(ind); 
        }
        void ValueRowInplace( const int & index, Eigen::MatrixXd &filledMatrix ) const { filledMatrix.row(index) = dat->row(index); }
        void ValueRowInplace( const int xIndex, const int yIndex, const int zIndex, Eigen::MatrixXd &filledMatrix ) const 
        { 
            const int ind = FlattenIndicies(xIndex, yIndex, zIndex);
            ValueRowInplace(ind,filledMatrix); 
        }
        void ValueScalarInplace( const int & index, T &filledBuffer ) const { assert(NData()==1); filledBuffer = dat->row(index); }
        void ValueScalarInplace( const int xIndex, const int yIndex, const int zIndex, T &filledBuffer ) const 
        { 
            const int ind = FlattenIndicies(xIndex, yIndex, zIndex);
            ValueScalarInplace(ind,filledBuffer); 
        }
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& GetMatrix() { return *dat; } 

        inline int NCells() const { return dat->rows(); }
        inline int NData() const { return dat->cols(); }  // for scalar its 1
        inline int FlattenIndicies(const int x, const int y, const int z) const { return x+y*owner_.cellCount_[0]+z*owner_.cellCount_[0]*owner_.cellCount_[1]; }

        /// IO methods
        void PrintData() const
        {
            std::cout << "grid data size: " << dat->rows() << "x" << dat->cols() << std::endl; 
            std::cout << *dat << std::endl; 
        }

        friend std::ostream& operator << ( std::ostream & os, const GridData& data )
        { 
            os << "grid data name: " << data.name << std::endl
               << "          size: " << data.dat->rows() << "x" << data.dat->cols();
            return os; 
        }
};

/* 
 * Rectangular grid base and handle interpolation
 */
template <class T> 
class Grid
{
    public: 
        enum INTERP_SCHEME { trilinear=0 };

    protected: 
        /// bounding box of the grid
        Eigen::Vector3d minBound_; 
        Eigen::Vector3d maxBound_; 

        /// cell count in each directions
        Eigen::Vector3i cellCount_; 

        /// data on the grid stored with a map. 
        std::unordered_map<std::string, GridData<T> > cellCenteredData_; 
        std::unordered_map<std::string, GridData<T> >       vertexData_; 

        /// interpolation scheme
        INTERP_SCHEME interpolationScheme_; 

    public:
        /// constructors 
        Grid() {} 
        Grid(const Eigen::Vector3d & minBound, const Eigen::Vector3d & maxBound, const Eigen::Vector3i cellCount) : 
            minBound_(minBound), maxBound_(maxBound), 
            cellCount_(cellCount), interpolationScheme_(trilinear) { }
        Grid(const Grid& grid) : 
            minBound_(grid.minBound_), maxBound_(grid.maxBound_), 
            cellCount_(grid.cellCount_), interpolationScheme_(grid.interpolationScheme_) { }


        /// destructors 
        ~Grid(){ } 

        /// data IO
        void InsertCellCenteredData( const std::string name, std::shared_ptr<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>> dataIn )
        {
            if ( dataIn->rows() != cellCount_[0]*cellCount_[1]*cellCount_[2] )
            {
                std::string errmsg  = "**ERROR** wrong cell-centered data size: " + std::to_string(dataIn->rows()) + "x" + std::to_string(dataIn->cols()); 
                            errmsg += ", because it is defined on grid of cell count " + std::to_string(cellCount_[0]) + "x"
                                                                                       + std::to_string(cellCount_[1]) + "x"
                                                                                       + std::to_string(cellCount_[2]); 
                throw std::runtime_error( errmsg );
            }

            GridData<T> gridData( name, dataIn, *this ); 
            std::pair<std::string,GridData<T>& > mapData (name, gridData);
            cellCenteredData_.insert( mapData );
        }

        void InsertVertexData( const std::string & name, std::shared_ptr<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>> dataIn )
        {
            if ( dataIn->rows() != (cellCount_[0])*(cellCount_[1])*(cellCount_[2]) )
            {
                std::string errmsg  = "**ERROR** wrong vertex data size: " + std::to_string(dataIn->rows()) + "x" + std::to_string(dataIn->cols()); 
                            errmsg += ", because it is defined on grid of cell count " + std::to_string(cellCount_[0]) + "x"
                                                                                       + std::to_string(cellCount_[1]) + "x"
                                                                                       + std::to_string(cellCount_[2]); 
                            throw std::runtime_error( errmsg );
            }

            GridData<T> gridData( name, dataIn, *this ); 
            std::pair<std::string,GridData<T>& > mapData (name, gridData);
            vertexData_.insert( mapData );
        }

        inline int FlattenIndicies(const int x, const int y, const int z) const { return x+y*cellCount_[0]+z*cellCount_[0]*cellCount_[1]; }
        inline void DeleteCellCenteredData( const std::string & name_key ) { cellCenteredData_.erase( name_key ); }
        inline void DeleteVertexData( const std::string & name_key ) { vertexData_.erase( name_key ); }
        inline void UpdateCellCenteredData( const std::string &nameKey, std::shared_ptr<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>> dataIn)
        {
            DeleteCellCenteredData(nameKey); 
            InsertCellCenteredData(nameKey,dataIn); 
        }
        inline GridData<T>& GetCellCenteredData( const std::string &name ) { return cellCenteredData_.at( name ); }
        inline GridData<T>&       GetVertexData( const std::string &name ) { return       vertexData_.at( name ); }
        inline Eigen::Vector3i GetCellCount() const { return cellCount_;}

        /// interpolation helpers 
        inline void SetInterpolationScheme( const INTERP_SCHEME & interpScheme ) { interpolationScheme_ = interpScheme; } 
        inline INTERP_SCHEME GetInterpolationScheme() { return interpolationScheme_; } 

        friend GridData<T>; 

        friend std::ostream& operator << ( std::ostream& os, const Grid<T> grid ) 
        {
            os << "Grid"; 
            return os; 
        }



};



/* 
 * Rectangular uniform grid 
 */
template <class T> 
class UniformGrid : public Grid<T>
{

    public: 
        // how many gradient components needed
        enum CELL_GRADIENT_COMPONENT
        { ALL=0, X=1, Y=2, Z=3, XY=12, YZ=23, XZ=13 };

    protected: 
        /// some length scale we can reason about
        Eigen::Vector3d dimension_; 
        Eigen::Vector3d dx_; 

    public: 

        /// constructors 

        UniformGrid() { }

        UniformGrid(const Eigen::Vector3d& minBound, const Eigen::Vector3d& maxBound, const Eigen::Vector3i& cellCount) : 
            Grid<T>(minBound, maxBound, cellCount),
            dimension_(this->maxBound_ - this->minBound_) { 
                dx_ << dimension_[0]/(double)this->cellCount_[0], 
                dimension_[1]/(double)this->cellCount_[1], 
                dimension_[2]/(double)this->cellCount_[2]; 
            }

        /// copy constructors
        UniformGrid(const UniformGrid<T>& grid) : Grid<T>(grid) {
            dimension_ = this->maxBound_ - this->minBound_;
            dx_ << dimension_[0]/(double)this->cellCount_[0], 
                dimension_[1]/(double)this->cellCount_[1], 
                dimension_[2]/(double)this->cellCount_[2]; 
        }

        /// return cell center position
        virtual Eigen::Vector3d GetCellCenterPosition( const int ii, const int jj, const int kk ) const 
        {
            return Eigen::Vector3d(this->minBound_[0]+((double)(ii)+0.5)*dx_[0], 
                    this->minBound_[1]+((double)(jj)+0.5)*dx_[1], 
                    this->minBound_[2]+((double)(kk)+0.5)*dx_[2]); 
        }

        /// interpolate the vertex data with given position
        Eigen::VectorXd InterpolateVertexData( const std::string& dataName, const Eigen::Vector3d& position ) 
        {
            const GridData<T>& data = this->GetVertexData( dataName ); 

            Eigen::VectorXd value; 

            switch (this->interpolationScheme_)
            {
                case this->trilinear:

                    value = TrilinearInterpolate<T>( this->minBound_, this->maxBound_, this->cellCount_, position, data );

                    break; 
            }

            return value; 
        }

        /// interpolate the cell-centered data with given position
        Eigen::VectorXd InterpolateCellCenteredData( const std::string& dataName, const Eigen::Vector3d& position )
        {

            const GridData<T>& data = this->GetCellCenteredData( dataName ); 

            Eigen::VectorXd value; 

            switch (this->interpolationScheme_)
            {
                case this->trilinear:

                    // define a few useful quantities for cell center
                    // interpolation
                    Eigen::Vector3d cmin = this->minBound_ + dx_/2.0; 
                    Eigen::Vector3d cmax = this->maxBound_ - dx_/2.0; 
                    Eigen::Vector3i Nc = this->cellCount_ - Eigen::Vector3i::Ones(3); 
                    value = TrilinearInterpolate<T>( cmin, cmax, Nc, position, data );

                    break; 
            }

            return value; 

        }

        // this is a naive implementation of the smoothing operation using
        // filter by N^3 iteration. 
        //
        // filter is passed in as an odd vector, the actual volumetric filter
        // will be constructed in the function. e.g.
        //
        // filter = [ 1 2 1 ] ---> volumertic filter
        // = [ 1 2 1 ] x [ 1 2 1 ] x [ 1 2 1 ]
        //
        // therefore this function only supports separable filters. 
        //
        // TODO
        // only works for scalar data for now
        //
        // TODO
        // if separable filters are used, the complexity can be brought down by
        // using smart caching, see e.g.
        // http://blogs.mathworks.com/steve/2006/10/04/separable-convolution/
        Eigen::MatrixXd CellCenteredSmoothing( const std::string &dataName, const Eigen::VectorXd &filter) 
        {
            const GridData<T> &data = this->GetCellCenteredData( dataName ); 
            const int dataDimension = data.NData(); 
            const int N_cells = data.NCells(); 
            assert( dataDimension == 1 ); 
            assert( filter.size() %2 == 1);

            const int kernelHalfWidth = (filter.size()-1)/2;

            const int Nx = this->cellCount_[0]; 
            const int Ny = this->cellCount_[1]; 
            const int Nz = this->cellCount_[2]; 

            Eigen::MatrixXd smoothedData(N_cells,dataDimension); 
            smoothedData.setZero(); 

            int count = 0; 

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
                        const int flattenedIndex = data.FlattenIndicies(ii,jj,kk); 

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

        // in case no buffer is wanted
        void CellCenteredScalarHessian(const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &hessianData)
        {
            std::vector<std::shared_ptr<Eigen::MatrixXd>> gradientBuffer; 
            CellCenteredScalarHessian(dataName, gradientBuffer, hessianData);
        }

        // compute the hessian of a scalar data pointed by dataName
        // Assume symmetry, the returned matrix will store in the order: 
        // xx, xy, xz, yy, yz, zz
        //
        // to minimize reallocation of memory, can pass in an optionally empty pre-allocated
        // memory block for gradient
        void CellCenteredScalarHessian(const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &gradientBuffer, std::vector<std::shared_ptr<Eigen::MatrixXd>> &hessianData) 
        {

            const GridData<T> &gridData = this->GetCellCenteredData(dataName); 
            const int dataDimension = gridData.NData(); 
            const int N_cells = gridData.NCells(); 
            assert( dataDimension==1 ); 

            const int hessianComponents = 6;


            //std::vector<std::shared_ptr<Eigen::MatrixXd>> gradientBuffer; 

            CellCenteredDataGradient( dataName, gradientBuffer, ALL); 

            std::shared_ptr<Eigen::MatrixXd> gx = gradientBuffer[0]; 
            std::shared_ptr<Eigen::MatrixXd> gy = gradientBuffer[1]; 
            std::shared_ptr<Eigen::MatrixXd> gz = gradientBuffer[2]; 

            this->InsertCellCenteredData(dataName+"gx", gx);
            this->InsertCellCenteredData(dataName+"gy", gy);
            this->InsertCellCenteredData(dataName+"gz", gz);

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

            this->DeleteCellCenteredData(dataName+"gx"); 
            this->DeleteCellCenteredData(dataName+"gy"); 
            this->DeleteCellCenteredData(dataName+"gz"); 

        }

        // take gradient to the data using finite-difference 
        // for example, for scalar data this will produce a vector3 at each
        // point. for vector3 data it will create a 3x3 matrix. 
        // 
        // for now only scalar data is tested. 
        //
        // 2nd finite-difference for 1st-order derivative 
        //   (f_i+1 - f_i-1) / 2h,              if 0<i<N-1
        //   (-3f_i + 4f_i+1 - f_i+2) / 2h,     if i=0
        //   (+3f_i - 4f_i-1 + f_i-2) / 2h,     if i=N-1
        //
        // returned flatted matrix consists on gradient components on the
        // columns and depends on component argument: 
        //
        // ALL: [gx, gy, gz] 
        // X:   [gx]
        // Y:   [gy]
        // Z:   [gz]
        // XY:  [gx,gy]
        // YZ:  [gy,gz]
        // XZ:  [gx,gz]
        //     
        void CellCenteredDataGradient( const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &gradientData , const CELL_GRADIENT_COMPONENT &component=ALL)
        {
            const GridData<T> & data = this->GetCellCenteredData( dataName ); 
            const int dataDimension =  data.NData(); 
            const int NCell = data.NCells(); 
            assert( dataDimension == 1 ); 

            const int Nx = this->cellCount_[0]; 
            const int Ny = this->cellCount_[1]; 
            const int Nz = this->cellCount_[2]; 

            const double dx = ( this->maxBound_[0] - this->minBound_[0] )/(double)Nx; 
            const double dy = ( this->maxBound_[1] - this->minBound_[1] )/(double)Ny; 
            const double dz = ( this->maxBound_[2] - this->minBound_[2] )/(double)Nz; 

            int N_gradientNeeded; 
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
                for (int ii=0; ii<N_gradientNeeded; ii++) 
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

        // please don't use this..
        std::vector<std::shared_ptr<Eigen::MatrixXd>> CellCenteredDataGradient( const std::string &dataName, const CELL_GRADIENT_COMPONENT &component=ALL ) 
        {
            std::vector<std::shared_ptr<Eigen::MatrixXd>> gradientData; 
            CellCenteredDataGradient(dataName, gradientData, component); 
            return gradientData; 
        }

        void PrintDataMeta() 
        {
            std::cout << "Print all data on the grid:\n"; 
            std::cout << "---------------------------\n"; 


            std::cout << "Cell-centered:\n";
            int count_cc=0, count_v=0; 
            std::cout << this->cellCenteredData_.size() << std::endl;
            //for (auto &p : this->cellCenteredData_) 
            //{
            //    std::cout << count_cc << " " << p.second << std::endl;
            //    count_cc ++; 
            //}

            //std::cout << "Vertex:\n";
            //for (auto &p : this->vertexData_) 
            //{
            //    std::cout << p.second << std::endl;
            //    count_v ++; 
            //}
        }

        void WriteVTKCellGrid(const std::string &vtkName)
        {
            const Eigen::Vector3i cellCount = this->GetCellCount(); 
            Eigen::MatrixXd gridPosition(cellCount[0]*cellCount[1]*cellCount[2],3); 

            for (int kk=0; kk<cellCount[2]; kk++) 
            {
                for (int jj=0; jj<cellCount[1]; jj++)
                {
                    for (int ii=0; ii<cellCount[0]; ii++)
                    {
                        const int index  = this->FlattenIndicies(ii,jj,kk); 
                        gridPosition.row(index) = this->GetCellCenterPosition(ii,jj,kk);
                    }
                }
            }

            VTKConverter::VTKStructureGridFromEigen( gridPosition, vtkName, VTKConverter::BINARY, cellCount );
        }

        void WriteVTKCellCentered(const std::string &vtkPrefix, const std::string &key, const std::string &dataName) 
        {
            const Eigen::Vector3i cellCount = this->GetCellCount(); 
            Eigen::MatrixXd gridPosition(cellCount[0]*cellCount[1]*cellCount[2],3); 
            GridData<double> &gridData = this->GetCellCenteredData(key);
            Eigen::MatrixXd &data = gridData.GetMatrix(); 

            const int N_data = data.cols(); 

            for (int kk=0; kk<cellCount[2]; kk++) 
            {
                for (int jj=0; jj<cellCount[1]; jj++)
                {
                    for (int ii=0; ii<cellCount[0]; ii++)
                    {
                        const int index  = gridData.FlattenIndicies(ii,jj,kk); 
                        gridPosition.row(index) = this->GetCellCenterPosition(ii,jj,kk);
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

        void WriteCSVCellCentered(const std::string &csvName, const std::string &dataName)
        {
            std::ofstream of(csvName.c_str());
            const int &Nx = this->cellCount_[0]; 
            const int &Ny = this->cellCount_[1]; 
            const int &Nz = this->cellCount_[2]; 
            const Eigen::Vector3d &maxBound = this->maxBound_;
            const Eigen::Vector3d &minBound = this->minBound_;
            const Eigen::Vector3d &division = this->division_;
            Eigen::MatrixXd &data = this->GetCellCenteredData(dataName); 

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


        friend std::ostream& operator << ( std::ostream& os, const UniformGrid<T> grid ) 
        { 

            os << std::setprecision(4) << std::fixed;
            os << "Uniform grid: \n"; 
            os << "------------- \n"; 
            os << "min: [" << grid.minBound_.transpose()  << "]\n";  
            os << "max: [" << grid.maxBound_.transpose()  << "]\n";  
            os << "dim: [" << grid.dimension_.transpose() << "]\n";
            os << "dx : [" << grid.dx_.transpose() << "]\n"; 
            os << "cell count : [" << grid.cellCount_.transpose() << "]\n";

            return os; 
        }

};



template <class T> 
class RotatedUniformGrid : public UniformGrid<T> 
{
    protected: 
        double yaw_;  // rotation about y-axis in radians 
        Eigen::Matrix3d rotation_; 

    public: 
        RotatedUniformGrid(const Eigen::Vector3d &minBound, const Eigen::Vector3d &maxBound, const Eigen::Vector3i &cellCount) : 
            UniformGrid<T>(minBound, maxBound, cellCount),
            yaw_(0) {
            }

        inline void   SetYawRadian(const double &yaw) { yaw_ = yaw; SetRotationMatrix();} 
        inline double GetYawRadian() { return yaw_; } 
        inline void   SetYawDegree(const double &yaw) { SetYawRadian(yaw/180.*M_PI); } 
        inline double GetYawDegree() { return yaw_*180./M_PI; } 

        void SetRotationMatrix()
        {
            Eigen::Matrix3d rotation; 
            rotation << cos(yaw_), 0,       sin(yaw_), 
                        0        , 1,       0, 
                       -sin(yaw_), 0,       cos(yaw_); 
            rotation_ = rotation; 
        }

        virtual Eigen::Vector3d GetCellCenterPosition( const int ii, const int jj, const int kk ) const 
        {

            Eigen::Vector3d position(this->minBound_[0]+((double)(ii)+0.5)*this->dx_[0], 
                                     this->minBound_[1]+((double)(jj)+0.5)*this->dx_[1],
                                     this->minBound_[2]+((double)(kk)+0.5)*this->dx_[2]); 

            return rotation_*position; 

        }



}; 

//-----------------------------------------------------------------------------------------------------------------------------------------//


/*
 * Helper for trilinear interpolation on the grid.
 *
 * trilinear interpolation following notation on wiki
 * (https://en.wikipedia.org/wiki/Trilinear_interpolation)
 *
 * this implementation does not depend on grid size
 */
    template <class T> 
Eigen::VectorXd TrilinearInterpolate( const Eigen::Vector3d& xmin, const Eigen::Vector3d& xmax, const Eigen::Vector3i& N, const Eigen::Vector3d& position, const GridData<T>& data ) 
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
    x = std::max<T>( std::min<T>( position[0], xmax[0] ), xmin[0] );
    y = std::max<T>( std::min<T>( position[1], xmax[1] ), xmin[1] );
    z = std::max<T>( std::min<T>( position[2], xmax[2] ), xmin[2] );
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

        T c00 = data.Value(ix0,iy0,iz0)(ii)*(1.0-xd) + data.Value(ix1,iy0,iz0)(ii)*xd; 
        T c10 = data.Value(ix0,iy1,iz0)(ii)*(1.0-xd) + data.Value(ix1,iy1,iz0)(ii)*xd; 
        T c01 = data.Value(ix0,iy0,iz1)(ii)*(1.0-xd) + data.Value(ix1,iy0,iz1)(ii)*xd; 
        T c11 = data.Value(ix0,iy1,iz1)(ii)*(1.0-xd) + data.Value(ix1,iy1,iz1)(ii)*xd; 

        T c0 = c00*(1.0 - yd) + c10*yd; 
        T c1 = c01*(1.0 - yd) + c11*yd; 

        value[ii] = c0*(1.0 - zd) + c1*zd; 
    }

    return value;
}


#endif 
