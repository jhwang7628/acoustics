#ifndef CUDA_2D_WAVE_SOLVER_H
#define CUDA_2D_WAVE_SOLVER_H

#include "../defines.h"

struct Cuda_Wave_2d_sim_data_t;

typedef Cuda_Wave_2d_sim_data_t* Cuda_Wave_2d_t;

Cuda_Wave_2d_t wave_sim_init(Number_t xmin, Number_t ymin,
							Number_t xmax, Number_t ymax,
							Number_t c, Number_t dt,
							int nx, int ny,
							Number_t (*init_function)(Number_t, Number_t, void *),
							void * ctx);

void wave_sim_free(Cuda_Wave_2d_t wave);

Number_t wave_sim_get_x(Cuda_Wave_2d_t wave, int j);
Number_t wave_sim_get_y(Cuda_Wave_2d_t wave, int i);

//Acessing element:
// u[j + i*(nx+2)]
// i -> y, j->x

Number_t * wave_sim_get_u(Cuda_Wave_2d_t wave, int offset);

void wave_sim_step(Cuda_Wave_2d_t wave);
void wave_sim_apply_boundary(Cuda_Wave_2d_t wave);

#endif