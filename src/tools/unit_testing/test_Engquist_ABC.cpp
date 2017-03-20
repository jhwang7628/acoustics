#include <iostream> 
#include <iomanip>
#include "Eigen/Dense" 

//##############################################################################
// Typedefs and global var
//##############################################################################
typedef double T; 
typedef Eigen::MatrixXd Matrix2D;
static const T AIR_DENSITY = 1.2; 
static const T SOUND_SPEED = 343.; 
static const int N = 21;
static const T CELL_SIZE = 0.01;
static const T STEP_SIZE = 1./59409.;

//##############################################################################
// Forward declaration
//##############################################################################
class Grid; 
class Solver; 

//##############################################################################
// Class Grid
//##############################################################################
class Grid 
{
private: 
    std::vector<Matrix2D> _data; 
    int _p; 
public: 
    Grid(const int &Nx, const int &Ny) 
        : _p(1) 
    {
        _data.resize(3); 
        for (int ii=0; ii<3; ++ii) 
        {
            _data[ii].setZero(Nx, Ny); 
        }
    }
    inline Matrix2D &GetData() 
    {
        return _data.at(_p); 
    }
    inline void SetData(const int &ind_x, const int &ind_y, const T &d)
    {
        _data.at(_p)(ind_x, ind_y) = d; 
    }
friend Solver;
}; 

//##############################################################################
// Class Solver
//##############################################################################
class Solver 
{
private: 
    Grid &_grid;
public: 
    Solver(Grid &grid)
        : _grid(grid) 
    {}
    void Step(); 
};

//##############################################################################
// Function main
//##############################################################################
int main() 
{
    Grid grid(N, N); 
    Solver solver(grid); 
    grid.SetData(10, 10, 1.);
    int c = 0;
    while (c<20)
    {
        std::cout << "---------- step " << c << " ----------" << std::endl;
        std::cout << std::setprecision(3) << std::fixed 
                  << grid.GetData() << std::endl;
        solver.Step();
        ++c; 
    }
}

//##############################################################################
// Function Step
//##############################################################################
void Solver::
Step() 
{
    auto &data_n = _grid._data[(_grid._p + 1) % 3]; 
    auto &data_c = _grid._data[(_grid._p    ) % 3]; 
    auto &data_p = _grid._data[(_grid._p + 2) % 3]; 
    const T lambda  = SOUND_SPEED * STEP_SIZE / CELL_SIZE; 
    const T lambda2 = pow(lambda, 2); 
    for (int jj=0; jj<N; ++jj) 
    {
        for (int ii=0; ii<N; ++ii)
        {
            const int jj_p = std::min(jj+1, N-1); 
            const int jj_n = std::max(jj-1, 0  ); 
            const int ii_p = std::min(ii+1, N-1); 
            const int ii_n = std::max(ii-1, 0  ); 
            data_n(ii,jj) = 2.0*data_c(ii,jj) - data_p(ii,jj)
                          + lambda2*(data_c(ii_p,jj  )+data_c(ii_n,jj  ) 
                                    +data_c(ii  ,jj_p)+data_c(ii  ,jj_n)
                                    -4.0*data_c(ii,jj));
        }
    }
    _grid._p = (_grid._p + 1) % 3; 
}

