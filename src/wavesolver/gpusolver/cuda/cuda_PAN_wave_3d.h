#ifndef Cuda_PAN_3D_WAVE_SOLVER_H
#define Cuda_PAN_3D_WAVE_SOLVER_H

#include "../defines.h" 
#include <boost/function.hpp>

struct Cuda_PAN_Wave_3d_sim_data_t;

typedef Cuda_PAN_Wave_3d_sim_data_t * Cuda_PAN_Wave_3d_t;

Cuda_PAN_Wave_3d_t wave_sim_init(Number_t xmin, Number_t ymin, Number_t zmin,
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
								 Number_t pulse);

void wave_sim_step(Cuda_PAN_Wave_3d_t wave);
Number_t * wave_listen(Cuda_PAN_Wave_3d_t wave, int field);

void wave_sim_get_divisions(const Cuda_PAN_Wave_3d_t wave, int * nx, int * ny, int * nz);
Number_t wave_sim_get_current_time(const Cuda_PAN_Wave_3d_t wave);
void wave_sim_get_bounds(const Cuda_PAN_Wave_3d_t wave,
						 Number_t * xmin, Number_t * xmax,
						 Number_t * ymin, Number_t * ymax,
						 Number_t * zmin, Number_t * zmax);

void wave_sim_free(Cuda_PAN_Wave_3d_t wave);

Number_t wave_sim_get_x(const Cuda_PAN_Wave_3d_t wave, int j);
Number_t wave_sim_get_y(const Cuda_PAN_Wave_3d_t wave, int i);
Number_t wave_sim_get_z(const Cuda_PAN_Wave_3d_t wave, int k);

Number_t * wave_sim_get_u(Cuda_PAN_Wave_3d_t wave);

#endif