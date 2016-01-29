#ifndef GRID_H
#define GRID_H 

#include <iostream>
#include <string>
#include <unordered_map>
#include <iomanip> 
#include <memory>

#include <Eigen/Dense> 

#ifndef REAL 
#define REAL double 
#endif

class Grid;

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
class GridData
{
    private: 
        std::string name;
        std::shared_ptr<Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>> dat; 

        /// the grid that owns this data
        Grid& owner_; 

    public: 

        GridData( const std::string& nameIn, std::shared_ptr<Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>> dataIn, Grid& owner ) : name(nameIn), dat(dataIn), owner_(owner) { } 

        /// getters 
        inline Eigen::Matrix<REAL, 1, Eigen::Dynamic> Value( const int & index ) const { return dat->row(index); }
        inline Eigen::Matrix<REAL, 1, Eigen::Dynamic> Value( const int xIndex, const int yIndex, const int zIndex ) const { return dat->row(FlattenIndicies(xIndex, yIndex, zIndex)); }
        inline void ValueRowInplace( const int & index, Eigen::MatrixXd &filledMatrix ) const { filledMatrix.row(index) = dat->row(index); }
        inline void ValueRowInplace( const int xIndex, const int yIndex, const int zIndex, Eigen::MatrixXd &filledMatrix ) const { ValueRowInplace(FlattenIndicies(xIndex, yIndex, zIndex),filledMatrix); }
        inline void ValueScalarInplace( const int & index, REAL &filledBuffer ) const { assert(NData()==1); filledBuffer = (*dat)(index,0); }
        inline void ValueScalarInplace( const int xIndex, const int yIndex, const int zIndex, REAL &filledBuffer ) const { ValueScalarInplace(FlattenIndicies(xIndex, yIndex, zIndex),filledBuffer); }
        inline Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>& GetMatrix() { return *dat; } 
        inline int NCells() const { return dat->rows(); }
        inline int NData() const { return dat->cols(); }  // for scalar its 1
        int FlattenIndicies(const int x, const int y, const int z) const;

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
        std::unordered_map<std::string, GridData> cellCenteredData_; 
        std::unordered_map<std::string, GridData>       vertexData_; 

        /// interpolation scheme
        INTERP_SCHEME interpolationScheme_; 

    public:
        /// constructors 
        Grid() {} 
        Grid(const Eigen::Vector3d & minBound, const Eigen::Vector3d & maxBound, const Eigen::Vector3i cellCount) : 
            minBound_(minBound), 
            maxBound_(maxBound), 
            cellCount_(cellCount), 
            interpolationScheme_(trilinear) { }
        Grid(const Grid& grid) : 
            minBound_(grid.minBound_), 
            maxBound_(grid.maxBound_), 
            cellCount_(grid.cellCount_), 
            interpolationScheme_(grid.interpolationScheme_) { }

        /// destructors 
        ~Grid(){ } 

        inline void FlattenIndicies(const int &x, const int &y, const int &z, int &ind) const { ind = x+y*cellCount_[0]+z*cellCount_[0]*cellCount_[1]; } 
        inline int FlattenIndicies(const int &x, const int &y, const int &z) const { int ind; FlattenIndicies(x,y,z,ind); return ind; }
        inline void UnflattenIndicies(const int &ind, int &x, int &y, int &z) const 
        {
            x = ind      % cellCount_[0]; 
            y = ((ind-x) / cellCount_[0]) % cellCount_[1];
            z = (ind-x-y)/ cellCount_[0]/cellCount_[1]; 
        }
        inline void UnflattenIndicies(const int &ind, Eigen::Vector3i &indicies) const { UnflattenIndicies(ind,indicies[0],indicies[1],indicies[2]); }

        inline void DeleteCellCenteredData( const std::string & name_key ) { cellCenteredData_.erase( name_key ); }
        inline void DeleteVertexData( const std::string & name_key ) { vertexData_.erase( name_key ); }
        inline void UpdateCellCenteredData( const std::string &nameKey, std::shared_ptr<Eigen::Matrix<REAL,Eigen::Dynamic,Eigen::Dynamic>> dataIn)
        {
            DeleteCellCenteredData(nameKey); 
            InsertCellCenteredData(nameKey,dataIn); 
        }
        inline GridData& GetCellCenteredData( const std::string &name ) { return cellCenteredData_.at( name ); }
        inline GridData&       GetVertexData( const std::string &name ) { return       vertexData_.at( name ); }
        inline Eigen::Vector3i GetCellCount() const { return cellCount_;}

        /// interpolation helpers 
        inline void SetInterpolationScheme( const INTERP_SCHEME & interpScheme ) { interpolationScheme_ = interpScheme; } 
        inline INTERP_SCHEME GetInterpolationScheme() const { return interpolationScheme_; } 
        inline int N_cells() const { return cellCount_[0]*cellCount_[1]*cellCount_[2]; } 

        virtual Eigen::Vector3d GetCellCenterPosition( const int &ii, const int &jj, const int &kk ) const = 0;
        /// data IO
        void InsertCellCenteredData( const std::string &name, std::shared_ptr<Eigen::Matrix<REAL,Eigen::Dynamic,Eigen::Dynamic>> dataIn ); 
        void InsertVertexData( const std::string &name, std::shared_ptr<Eigen::Matrix<REAL,Eigen::Dynamic,Eigen::Dynamic>> dataIn ); 

        friend GridData; 
        friend std::ostream& operator << ( std::ostream& os, const Grid &grid ) 
        {
            os << "Grid"; 
            return os; 
        }
};



/* 
 * Rectangular uniform grid 
 */
class UniformGrid : public Grid
{

    public: 
        // how many gradient components needed
        enum CELL_GRADIENT_COMPONENT { ALL=0, X=1, Y=2, Z=3, XY=12, YZ=23, XZ=13 };

    protected: 
        /// some length scale we can reason about
        Eigen::Vector3d dx_; 

    public: 

        /// constructors 
        UniformGrid() { }
        UniformGrid(const Eigen::Vector3d& minBound, const Eigen::Vector3d& maxBound, const Eigen::Vector3i& cellCount); 
        /// copy constructors
        UniformGrid(const UniformGrid& grid);

        inline void RecomputeCachedField() { dx_ = (maxBound_ - minBound_).cwiseQuotient(cellCount_.cast<double>()); } 
        virtual Eigen::Vector3d GetCellCenterPosition( const int &ii, const int &jj, const int &kk ) const; 
        virtual inline void GetCellCenterPosition(const int &index, Eigen::Vector3d &position) const {int ii,jj,kk; UnflattenIndicies(index,ii,jj,kk); GetCellCenterPosition(ii,jj,kk,position[0],position[1],position[2]);} 
        virtual inline void GetCellCenterPosition(const int &ii, const int &jj, const int &kk, Eigen::Vector3d &position) const {GetCellCenterPosition(ii,jj,kk,position[0],position[1],position[2]);} 
        virtual void GetCellCenterPosition( const int &ii, const int &jj, const int &kk, double &x, double &y, double &z ) const; 
        virtual void GetAllCellCenterPosition(Eigen::MatrixXd &gridPosition) const; 
        static void GetAllCellCenterPosition(const Eigen::Vector3d &minBound, const Eigen::Vector3d &maxBound, const Eigen::Vector3i &cellCount, Eigen::MatrixXd &gridPosition); 

        /// interpolate the vertex data with given position
        Eigen::VectorXd InterpolateVertexData( const std::string& dataName, const Eigen::Vector3d& position ); 

        /// interpolate the cell-centered data with given position
        Eigen::VectorXd InterpolateCellCenteredData( const std::string& dataName, const Eigen::Vector3d& position ); 

        // *filter is passed in as an odd vector, the actual volumetric filter
        //  will be constructed in the function. e.g.
        //
        //  filter = [ 1 2 1 ] ---> volumertic filter
        //  = [ 1 2 1 ] x [ 1 2 1 ] x [ 1 2 1 ]
        //
        //  therefore this function only supports separable filters. 
        // 
        // *only where the mask is true will be filtered. 
        //
        // TODO
        // only works for scalar data for now
        Eigen::MatrixXd CellCenteredSmoothing( const std::string &dataName, const Eigen::VectorXd &filter, const std::vector<bool> &mask);

        // in case no buffer is wanted
        virtual void CellCenteredScalarHessian(const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &hessianData); 

        // compute the hessian of a scalar data pointed by dataName
        // Assume symmetry, the returned matrix will store in the order: 
        // xx, xy, xz, yy, yz, zz
        //
        // to minimize reallocation of memory, can pass in an optionally empty pre-allocated
        // memory block for gradient
        virtual void CellCenteredScalarHessian(const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &gradientBuffer, std::vector<std::shared_ptr<Eigen::MatrixXd>> &hessianData);

        // take gradient to the data using finite-difference 
        // for example, for scalar data this will produce a vector3 at each
        // point. for vector3 data it will create a 3x3 matrix. 
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
        // TODO
        // for now only scalar data is tested. 
        virtual void CellCenteredDataGradient( const std::string &dataName, std::vector<std::shared_ptr<Eigen::MatrixXd>> &gradientData , const CELL_GRADIENT_COMPONENT &component=ALL); 

        void PrintDataMeta();

        void WriteVTKCellGrid(const std::string &vtkName); 
        void WriteVTKCellCentered(const std::string &vtkPrefix, const std::string &key, const std::string &dataName); 

        // static method to write vtk directly from eigen data
        static void WriteVTKCellCenteredFromEigen(std::shared_ptr<Eigen::MatrixXd> ptrData, std::shared_ptr<Eigen::MatrixXd> ptrGridPosition, const std::string &vtkPrefix, const std::string &dataName); 
        static void WriteVTKCellCenteredFromEigen(std::shared_ptr<Eigen::MatrixXd> ptrData, const Eigen::Vector3d &minBound, const Eigen::Vector3d &maxBound, const Eigen::Vector3i &cellCount, const std::string &vtkPrefix, const std::string &dataName); 

        // a simple wrapper that will copy the incoming data to write to vtk.
        // it will clean up after itself.  
        static void WriteVTKCellCenteredFromEigen(const Eigen::MatrixXd &data, const Eigen::Vector3d &minBound, const Eigen::Vector3d &maxBound, const Eigen::Vector3i &cellCount, const std::string &vtkPrefix, const std::string &dataName); 

        void WriteCSVCellCentered(const std::string &csvName, const std::string &dataName); 

        friend std::ostream& operator << ( std::ostream& os, const UniformGrid grid ) 
        { 

            os << std::setprecision(12) << std::fixed;
            os << "Uniform grid: \n"; 
            os << "------------- \n"; 
            os << "min: [" << grid.minBound_.transpose()  << "]\n";  
            os << "max: [" << grid.maxBound_.transpose()  << "]\n";  
            os << "dx : [" << grid.dx_.transpose() << "]\n"; 
            os << "cell count : [" << grid.cellCount_.transpose() << "]\n";

            return os; 
        }

};



class RotatedUniformGrid : public UniformGrid 
{
    protected: 
        double yaw_;  // rotation about y-axis in radians 
        Eigen::Matrix3d rotation_; 

    public: 
        RotatedUniformGrid(const Eigen::Vector3d &minBound, const Eigen::Vector3d &maxBound, const Eigen::Vector3i &cellCount) : 
            UniformGrid(minBound, maxBound, cellCount),
            yaw_(0) { }

        inline void   SetYawRadian(const double &yaw) { yaw_ = yaw; SetRotationMatrix();} 
        inline double GetYawRadian() { return yaw_; } 
        inline void   SetYawDegree(const double &yaw) { SetYawRadian(yaw/180.*M_PI); } 
        inline double GetYawDegree() { return yaw_*180./M_PI; } 

        void SetRotationMatrix(); 

        virtual Eigen::Vector3d GetCellCenterPosition( const int ii, const int jj, const int kk ) const; 
}; 

//-----------------------------------------------------------------------------------------------------------------------------------------//


// helper
Eigen::VectorXd TrilinearInterpolate( const Eigen::Vector3d& xmin, const Eigen::Vector3d& xmax, const Eigen::Vector3i& N, const Eigen::Vector3d& position, const GridData& data );


#endif 
