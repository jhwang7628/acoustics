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
static const T AIR_DENSITY = 1.2; 
static const T SOUND_SPEED = 343.; 
static const int N = 100;
static const T CELL_SIZE = 0.005;
// derived
static const T STEP_SIZE = CELL_SIZE / (sqrt(3.)*SOUND_SPEED);
static const T BOX_SIZE = CELL_SIZE*(T)N;

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
    Vector2 _minBound; 
    Vector2 _maxBound; 
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
        const T low_corner = -BOX_SIZE/2.;
        const T top_corner =  BOX_SIZE/2.;
        _minBound << low_corner, low_corner;
        _maxBound << top_corner, top_corner;
    }
    inline Vector2 CellPosition(const int &x, const int &y)
    {
        Vector2 pos = _minBound;
        pos.array() += CELL_SIZE/2.0;
        pos[0] += (T)x * CELL_SIZE;
        pos[1] += (T)y * CELL_SIZE;
        return pos;
    }
    inline Matrix2D &GetData() 
    {
        return _data.at(_p); 
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
int main(int argc, char **argv) 
{
    // parse
    const int N_steps = (argc==1 ? 200 : atoi(argv[1]));
    std::cout << STEP_SIZE << std::endl;
    // run
    Grid grid(N, N); 
    Solver solver(grid); 
    grid.InitializeGaussian(0.02);
    char filename[512];
    int c = 0;
    while (c<N_steps)
    {
        std::cout << "step " << c << "" << std::endl;
        snprintf(filename, 512, "data/%.5d.dat", c); 
        grid.SaveData(filename);
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

//##############################################################################
// Function SaveData
//##############################################################################
void Grid::
SaveData(const char *filename)
{
    std::ofstream stream(filename, std::ios::out|std::ios::binary); 
    stream.write((char*)(&N), sizeof(int)); // rows
    stream.write((char*)(&N), sizeof(int)); // cols
    stream.write((char*)(_data.at(_p).data()), sizeof(T)*N*N); //data
    stream.close(); 
}

//##############################################################################
// Function InitializeGaussian
//##############################################################################
void Grid::
InitializeGaussian(const T &r)
{
    auto Gaussian = [=](const Vector2 &pos)
    {
        return exp(-pow(pos[0]/r, 2) - pow(pos[1]/r, 2)); 
    };
    for (int jj=0; jj<N; ++jj)
    {
        for (int ii=0; ii<N; ++ii)
        {
            const Vector2 pos = CellPosition(ii, jj); 
            SetData(ii, jj, Gaussian(pos)); 
        }
    }
}
