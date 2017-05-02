#include <iostream> 
#include <fstream> 
#include <iomanip>
#include <string>
#include "Eigen/Dense" 

//##############################################################################
// Typedefs and global var
//##############################################################################
typedef double T; 
typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Matrix2D;
typedef Eigen::Matrix<T, 2, 1> Vector2;
typedef Eigen::Matrix<int, 2, 1> Vector2i;
static const T AIR_DENSITY = 1.2; 
static const T SOUND_SPEED = 343.; 
static const int N = 150;
static const T CELL_SIZE = 0.0025;
// derived
static const T STEP_SIZE = CELL_SIZE / (sqrt(3.)*SOUND_SPEED);
static const T BOX_SIZE = CELL_SIZE*(T)N;

//##############################################################################
// Helper functions
//##############################################################################
template <typename _T> 
_T clamp(const _T &low, const _T &high, const _T &val)
{
    return std::max(low, std::min(high, val)); 
}

//##############################################################################
// Forward declaration
//##############################################################################
class Grid; 
class WorldGrid; 
class ComputeGrid; 
class Solver; 
class DivergenceSource; 

//##############################################################################
// Class DivergenceSource
//##############################################################################
struct DivergenceSource
{
private: 
    T _radius_space; 
    T _radius_time; 
public: 
    Vector2 center; 
    T radius; 
    T start = 0.;
    inline void SetRadius(const T &r)
    {
        radius = r; 
        _radius_time = r; 
        _radius_space = SOUND_SPEED*_radius_time; 
        start = STEP_SIZE*(T)50; 
    }
    inline bool Inside(const Vector2 &evalpos, const T &t)
    {
        return ((evalpos - center).norm() < _radius_space && abs(t-start)<_radius_time); 
    }
    inline T Evaluate(const Vector2 &evalpos, const T &t)
    {
        T value = (evalpos-center).squaredNorm(); 
        value = -value / (2.0*pow(_radius_space,2)); 
        value = exp(value); 
        value *= -exp(-pow((t-start)/_radius_time,2)/2.0); 
        return value;
    }
};

//##############################################################################
// Class Grid
//##############################################################################
class Grid 
{
protected: 
    enum Type{COMPUTE, WORLD} _type;
    std::vector<Matrix2D> _data; 
    Vector2 _minBound; 
    Vector2 _maxBound; 
    int _p; 
    DivergenceSource *_source; 
public: 
    Grid(const int &Nx, const int &Ny) 
        : _p(1), _source(nullptr)
    {
        const T low_corner = -BOX_SIZE/2.;
        const T top_corner =  BOX_SIZE/2.;
        _minBound << low_corner, low_corner;
        _maxBound << top_corner, top_corner;
        _data.resize(3); 
        for (int ii=0; ii<3; ++ii) 
        {
            _data[ii].setZero(Nx, Ny); 
        }
    }
    ~Grid()
    {
        if (_source) delete _source; 
    }
    inline Vector2i Dimension() const
    {
        return {(int)_data[_p].rows(), (int)_data[_p].cols()}; 
    }
    inline Vector2 CellPosition(const int &x, const int &y) const
    {
        Vector2 pos = _minBound;
        pos.array() += CELL_SIZE/2.0;
        pos[0] += (T)x * CELL_SIZE;
        pos[1] += (T)y * CELL_SIZE;
        return pos;
    }
    inline Vector2i CellIndex(const Vector2 &pos) const
    {
        Vector2i ind; 
        ind[0] = clamp(0, (int)_data[_p].rows()-1, 
                       (int)((pos[0] - _minBound[0])/CELL_SIZE)); 
        ind[1] = clamp(0, (int)_data[_p].cols()-1, 
                       (int)((pos[1] - _minBound[1])/CELL_SIZE)); 
        return ind; 
    }
    inline std::vector<Matrix2D> &GetAllData()
    {
        return _data; 
    }
    inline Matrix2D &GetData(const int &p=-1) 
    {
        return (p<0 ? _data.at(_p) : _data.at(p)); 
    }
    inline void SetData(const int &ind_x, const int &ind_y, const T &d)
    {
        _data.at(_p)(ind_x, ind_y) = d; 
    }
    void SaveData(const char *filename); 
    void InitializeGaussian(const T &radius);
friend Solver;
}; 

//##############################################################################
// Class WorldGrid
//##############################################################################
class WorldGrid : public Grid
{
private: 
    static constexpr T NO_VAL = std::numeric_limits<T>::max();
    ComputeGrid *_computeGrid; 
public: 
    WorldGrid(const int &Nx, const int &Ny) 
        : Grid(Nx, Ny), _computeGrid(nullptr)
    {
        _type = WORLD; 
        _p = 0;
        ResetData();
    }
    inline void ResetData()
    {
        for (int ii=0; ii<3; ++ii)
        {
            _data[ii].setZero();
            //_data[ii].setOnes();
            //_data[ii].array()*=NO_VAL; 
        }
    }
    inline void SetComputeGrid(ComputeGrid *grid){_computeGrid = grid;}
    void MoveComputeGrid(const int &direction, const T &amount); 
};

//##############################################################################
// Class ComputeGrid
//##############################################################################
class ComputeGrid : public Grid
{
public: 
    ComputeGrid(const int &Nx, const int &Ny)
        : Grid(Nx, Ny)
    {
        _type = COMPUTE; 
    }
friend WorldGrid; 
};

//##############################################################################
// Class Solver
//##############################################################################
class Solver 
{
private: 
    Grid &_grid;
    T _time;
public: 
    Solver(Grid &grid)
        : _grid(grid), _time(0)
    {}
    void Step(); 
};

//##############################################################################
// Function main
//##############################################################################
int main(int argc, char **argv) 
{
    // parse
    const int N_steps = (argc==1 ? 200 : atoi(argv[1]));
    // run
    Grid *cgrid = new ComputeGrid(N+1, N); 
    WorldGrid *wgrid = new WorldGrid(2*N+1, N); 
    Solver solver(*cgrid); 
    cgrid->InitializeGaussian(STEP_SIZE*4.);
    wgrid->SetComputeGrid((ComputeGrid*)cgrid); 
    char filename[512];
    int c = 0;
    while (c<N_steps)
    {
        std::cout << "step " << c << "" << std::endl;
        snprintf(filename, 512, "data/%.5d.dat", c); 
        wgrid->SaveData(filename);
        //wgrid->MoveComputeGrid(0, CELL_SIZE/10.);
        wgrid->MoveComputeGrid(0, 0.);
        solver.Step();
        ++c; 
    }
    std::cout << "finish simulation with N=" << N << "\n"; 
    std::cout << "data stored in : " << filename << "\n"; 
}

//##############################################################################
// Function Step
//##############################################################################
void Solver::
Step() 
{
    auto &data_n = _grid._data[(_grid._p + 1) % 3]; //next
    auto &data_c = _grid._data[(_grid._p    ) % 3]; //current
    auto &data_p = _grid._data[(_grid._p + 2) % 3]; //past
    const T lambda  = SOUND_SPEED * STEP_SIZE / CELL_SIZE; 
    const T lambda2 = pow(lambda, 2); 
    for (int jj=0; jj<N; ++jj) 
    {
        for (int ii=0; ii<N; ++ii)
        {
            const Vector2 position = _grid.CellPosition(ii, jj); 
            if (ii==N-1)// || ii==1)
            {   
                const int jj_p = std::min(jj+1, N-1); 
                const int jj_n = std::max(jj-1, 0  ); 
                const int ii_p = (ii==N-1 ? ii+1 : ii-1); 
                const int ii_n = (ii==N-1 ? ii-1 : ii+1); 
                // upwinding
                //data_n(ii,jj) = data_c(ii,jj) - lambda*(data_c(ii_n,jj));  
                // Bilbao
                //data_n(ii,jj) = lambda2/(1.+lambda)*(
                //        (2./lambda2 - 4.)*data_c(ii,jj)
                //        + data_c(ii,jj_p) + data_c(ii,jj_n)
                //        +2.*data_c(ii_n,jj) - 1./lambda2*data_p(ii,jj));
                // 1st order Engquist (using 1st order space discretization)
                //std::cout << "mine\n";
                //data_n(ii,jj) = lambda2/(1.+2.*lambda)*(
                //        (2./lambda2 + 2./lambda - 4.)*data_c(ii,jj)
                //        + data_c(ii,jj_p) + data_c(ii,jj_n)
                //        +2.*data_c(ii_n,jj) - 1./lambda2*data_p(ii,jj));
                // 1st order Engquist (using 2nd order space discretization)
                data_n(ii,jj) = lambda2/(1.+ lambda)*(
                        (2./lambda2 - 4.)*data_c(ii,jj)
                        + data_c(ii,jj_p) + data_c(ii,jj_n)
                        +2.*data_c(ii_n,jj) +(lambda - 1.0)/lambda2*data_p(ii,jj));
                data_n(ii_p,jj) = 2./lambda*(data_c(ii,jj)-data_n(ii,jj))
                                + data_c(ii_n,jj); 
                continue; 
            }
            const int jj_p = std::min(jj+1, N-1); 
            const int jj_n = std::max(jj-1, 0  ); 
            const int ii_p = ii+1; //std::min(ii+1, N-1); 
            const int ii_n = std::max(ii-1, 0  ); 
            data_n(ii,jj) = 2.0*data_c(ii,jj) - data_p(ii,jj)
                          + lambda2*(data_c(ii_p,jj  )+data_c(ii_n,jj  ) 
                                    +data_c(ii  ,jj_p)+data_c(ii  ,jj_n)
                                    -4.0*data_c(ii,jj));
            if (_grid._source && _grid._source->Inside(position, _time)) 
            {
                data_n(ii,jj) += -AIR_DENSITY*pow(SOUND_SPEED,2)*STEP_SIZE*
                                 _grid._source->Evaluate(position, _time);
            }
        }
    }
    // second order Engquist
    //for (int jj=0; jj<N; ++jj)
    //{
    //    const int jj_p = std::min(jj+1, N-1); 
    //    const int jj_n = std::max(jj-1, 0  ); 
    //    data_n(N, jj) =     (data_n(N-1,jj)-data_p(N-1,jj)+    data_p(N  ,jj))/(2.*CELL_SIZE*STEP_SIZE)
    //                  - 0.5*(               data_p(N  ,jj)-2.0*data_c(N  ,jj))/(   STEP_SIZE*STEP_SIZE)
    //                  - 0.5*(data_n(N-1,jj)+data_p(N-1,jj)-2.0*data_c(N-1,jj))/(   STEP_SIZE*STEP_SIZE)
    //                  + 0.25*(data_p(N  ,jj_p)+data_p(N  ,jj_n)-2.0*data_p(N  ,jj))/(CELL_SIZE*CELL_SIZE)
    //                  + 0.25*(data_p(N-1,jj_p)+data_p(N-1,jj_n)-2.0*data_p(N-1,jj))/(CELL_SIZE*CELL_SIZE); 
    //    data_n(N,jj) *= 2.*CELL_SIZE*STEP_SIZE*STEP_SIZE / (STEP_SIZE + CELL_SIZE); 
    //}
    //std::cout << "---------------" << std::endl;
    //std::cout << data_n.row(N) << std::endl;
    //std::cout << data_n.row(N-1) << std::endl;
    //std::cout << data_n.row(N-2) << std::endl;
    //std::cout << data_n.row(N-3) << std::endl;
    _grid._p = (_grid._p + 1) % 3; 
    _time += STEP_SIZE; 
}

//##############################################################################
// Function SaveData
//##############################################################################
void Grid::
SaveData(const char *filename)
{
    const int rows = _data[_p].rows(); 
    const int cols = _data[_p].cols(); 
    std::ofstream stream(filename, std::ios::out|std::ios::binary); 
    stream.write((char*)(&rows), sizeof(int)); // rows
    stream.write((char*)(&cols), sizeof(int)); // cols
    stream.write((char*)(_data.at(_p).data()), sizeof(T)*rows*cols); //data
    stream.close(); 
}

//##############################################################################
// Function InitializeGaussian
//##############################################################################
void Grid::
InitializeGaussian(const T &radius_time)
{
    //// divergence source
    _source = new DivergenceSource(); 
    _source->center.setZero(); 
    _source->SetRadius(radius_time); 
    // simply initialize gaussian
    //const T r = SOUND_SPEED*radius_time;
    //auto &data_n = _data[(_p + 1) % 3]; //next
    //auto &data_c = _data[(_p    ) % 3]; //current
    //auto Gaussian = [=](const Vector2 &pos)
    //{
    //    return exp(-pow(pos[0]/r, 2) - pow(pos[1]/r, 2)); 
    //};
    //for (int jj=0; jj<N; ++jj)
    //{
    //    for (int ii=0; ii<N; ++ii)
    //    {
    //        const Vector2 pos = CellPosition(ii, jj); 
    //        data_c(ii,jj) = Gaussian(pos); 
    //    }
    //}
}

//##############################################################################
// Function MoveComputeGrid
//##############################################################################
void WorldGrid::
MoveComputeGrid(const int &direction, const T &amount)
{
    if (!_computeGrid)
        return; 
    const Vector2i dim = _computeGrid->Dimension();
    ResetData(); 
    const Vector2i ind0 = CellIndex(_computeGrid->_minBound); 
    // move computegrid bounds
    (_computeGrid->_minBound)[direction] += amount; 
    (_computeGrid->_maxBound)[direction] += amount; 
    const Vector2i ind1 = CellIndex(_computeGrid->_minBound); 
    // transfer computegrid data to worldgrid and back
    for (int i=0; i<3; ++i)
    {
        _data[i].block(ind0[0], ind0[1], dim[0], dim[1]) = _computeGrid->GetData(i);
        _computeGrid->GetData(i) = _data[i].block(ind1[0], ind1[1], dim[0], dim[1]); 
    }
    _p = _computeGrid->_p;
}
