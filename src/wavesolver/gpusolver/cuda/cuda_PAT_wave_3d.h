#ifndef Cuda_PAT_3D_WAVE_SOLVER_H
#define Cuda_PAT_3D_WAVE_SOLVER_H

#include "../defines.h" 
#include <boost/function.hpp>

struct Cuda_PAT_Wave_3d_sim_data_t;

typedef Cuda_PAT_Wave_3d_sim_data_t * Cuda_PAT_Wave_3d_t;

Cuda_PAT_Wave_3d_t wave_sim_init(Number_t xmin, Number_t ymin, Number_t zmin,
								 Number_t xmax, Number_t ymax, Number_t zmax,
								 Number_t c, Number_t dt,
								 Number_t cellsize,
								 int listening_count,
								 Number_t * listening_positions,
								 const Wave_InitialCondition3D & initial,
								 const Wave_BoundaryEvaluator3D & boundary,
								 Number_t xcenter, Number_t ycenter, Number_t zcenter,
								 const Wave_GradientEvaluator3D & gradient,
								 Number_t pml_width,
								 Number_t pml_strength,
								 Number_t frequency,
								 int num_multipole_coef,
								 Number_t multipole_radius);

void wave_sim_step(Cuda_PAT_Wave_3d_t wave);
Number_t * wave_listen(Cuda_PAT_Wave_3d_t wave);
Number_t * wave_compute_multipole(Cuda_PAT_Wave_3d_t wave, Number_t radius);
void wave_test_multipole(Cuda_PAT_Wave_3d_t wave);

void wave_sim_get_divisions(const Cuda_PAT_Wave_3d_t wave, int * nx, int * ny, int * nz);
Number_t wave_sim_get_current_time(const Cuda_PAT_Wave_3d_t wave);
void wave_sim_get_bounds(const Cuda_PAT_Wave_3d_t wave,
						 Number_t * xmin, Number_t * xmax,
						 Number_t * ymin, Number_t * ymax,
						 Number_t * zmin, Number_t * zmax);

void wave_sim_free(Cuda_PAT_Wave_3d_t wave);

Number_t wave_sim_get_x(const Cuda_PAT_Wave_3d_t wave, int j);
Number_t wave_sim_get_y(const Cuda_PAT_Wave_3d_t wave, int i);
Number_t wave_sim_get_z(const Cuda_PAT_Wave_3d_t wave, int k);

bool checkIsBulk(const Cuda_PAT_Wave_3d_t wave, int i, int j, int k);
void wave_gradientAt(const Cuda_PAT_Wave_3d_t wave, int i, int j, int k, Number_t * x, Number_t * y, Number_t * z);

void wave_estimate_ijk(const Cuda_PAT_Wave_3d_t wave, Number_t x, Number_t y, Number_t z, int * i, int * j, int * k);
void wave_estimate_with_multipole(Cuda_PAT_Wave_3d_t wave, Number_t x, Number_t y, Number_t z, Number_t * amplitude, Number_t * phase);

Number_t * wave_sim_get_u(Cuda_PAT_Wave_3d_t wave);
Number_t * wave_sim_get_amplitudes(Cuda_PAT_Wave_3d_t wave);
Number_t * wave_sim_get_phases(Cuda_PAT_Wave_3d_t wave);

#endif
