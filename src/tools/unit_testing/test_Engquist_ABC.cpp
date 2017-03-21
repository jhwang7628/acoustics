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
static const int N = 150;
static const T CELL_SIZE = 0.005;
// derived
static const T STEP_SIZE = CELL_SIZE / (sqrt(3.)*SOUND_SPEED);
static const T BOX_SIZE = CELL_SIZE*(T)N;

//##############################################################################
// Forward declaration
//##############################################################################
class Grid; 
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
private: 
    std::vector<Matrix2D> _data; 
    Vector2 _minBound; 
    Vector2 _maxBound; 
    int _p; 
    DivergenceSource *_source; 
public: 
    Grid(const int &Nx, const int &Ny) 
        : _p(1), _source(nullptr)
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
    ~Grid()
    {
        if (_source) delete _source; 
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
    Grid grid(N+1, N); 
    Solver solver(grid); 
    grid.InitializeGaussian(STEP_SIZE*10.);
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
            if (ii==N-1 && jj!=0 && jj!=N-1)
            {   
                // first-order upwinding
                //data_n(ii,jj) = data_c(ii,jj) - lambda*(data_c(ii-1,jj));  
                // second-order Bilbao
                //data_n(ii,jj) = lambda2/(1.+lambda)*(
                //        (2./lambda2 - 4.)*data_c(ii,jj)
                //        + data_c(ii,jj+1) + data_c(ii,jj-1)
                //        +2.*data_c(ii-1,jj) - 1./lambda2*data_p(ii,jj));
                // second-order mine
                data_n(ii,jj) = lambda2/(1.+2.*lambda)*(
                        (2./lambda2 + 2./lambda - 4.)*data_c(ii,jj)
                        + data_c(ii,jj+1) + data_c(ii,jj-1)
                        +2.*data_c(ii-1,jj) - 1./lambda2*data_p(ii,jj));

                data_n(ii+1,jj) = 2./lambda*(data_c(ii,jj)-data_n(ii,jj))
                                + data_c(ii-1,jj); 

                continue; 
            }
            const int jj_p = std::min(jj+1, N-1); 
            const int jj_n = std::max(jj-1, 0  ); 
            const int ii_p = std::min(ii+1, N-1); 
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
    //for (int jj=1; jj<N-1; ++jj)
    //{
    //    for (int ii=1; ii<N-1; ++ii)
    //    {
    //        data_n(ii,jj) = data_c(ii,jj) + pow(SOUND_SPEED*STEP_SIZE/CELL_SIZE,2)
    //                                       *(data_c(ii+1,jj  )-data_c(ii,jj)+data_c(ii-1,jj  )
    //                                        +data_c(ii  ,jj+1)-data_c(ii,jj)+data_c(ii  ,jj-1)); 
    //    }
    //}
    //_p = (_p+1)%3;
}
