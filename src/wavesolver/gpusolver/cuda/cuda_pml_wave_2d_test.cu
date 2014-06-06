#include "cuda_pml_wave_2d.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cstdio>

#include "cuda_helper.h"

#include "cuda_pml_wave_2d_kernel_test.cu"

struct Cuda_PML_Wave_2d_sim_data_t {
	Number_t xmin, ymin;
	Number_t xmax, ymax;

	Number_t dt;
	Number_t dx;
	Number_t dy;

	Number_t t;

	Number_t c;

	int nx, ny;

	Number_t pml_width;
	Number_t pml_strength;
	Number_t density;

	Number_t * ubuf;
	Number_t * ubuf_d;
	bool * isBulk;
	bool * isBulk_d;
	Number_t * gradient;
	Number_t * gradient_d;
	Number_t radius;
};

Number_t * wave_sim_get_u(Cuda_PML_Wave_2d_t wave){
	return wave->ubuf;
}

Number_t gradient(Number_t x, Number_t y, int dim){
	Number_t dist = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
	if(dist < 0.33){
		if(dim == 0){
			return -(x-0.5)/dist;
		} else{
			return -(y-0.5)/dist;
		}
	} else{
		if(dim == 0){
			return (x-0.5)/dist;
		} else{
			return (y-0.5)/dist;
		}
	}
}

Cuda_PML_Wave_2d_t wave_sim_init(Number_t xmin, Number_t ymin,
						Number_t xmax, Number_t ymax,
						Number_t c, Number_t dt, 
						int nx, int ny,
						Number_t (*init_function)(Number_t, Number_t, void *),
						void * ctx,
						Number_t pml_width,
						Number_t pml_strength){
	Cuda_PML_Wave_2d_t wave = (Cuda_PML_Wave_2d_t) malloc(sizeof(Cuda_PML_Wave_2d_sim_data_t));
	
	wave->xmin = xmin;
	wave->ymin = ymin;
	wave->xmax = xmax;
	wave->ymax = ymax;
	wave->radius = 0.35;
	Number_t radius = wave->radius;
	wave->c = c;

	wave->dt = dt;

	wave->t = 0;
	wave->density = 1;

	wave->nx = nx;
	wave->ny = ny;

	wave->dx = (xmax-xmin)/nx;
	wave->dy = (ymax-ymin)/ny;

	wave->pml_width = pml_width;
	wave->pml_strength = pml_strength;
	wave->ubuf = NULL;
	cudaCheckError(cudaMallocHost((void**)&wave->ubuf, 3*(wave->nx)*(wave->ny)*sizeof(Number_t)));
	wave->isBulk = NULL;
	cudaCheckError(cudaMallocHost((void**)&wave->isBulk, (wave->nx)*(wave->ny)*sizeof(bool)));
	wave->gradient = NULL;
	cudaCheckError(cudaMallocHost((void**)&wave->gradient, 2*(wave->nx)*(wave->ny)*sizeof(Number_t)));
	
	wave->ubuf_d = NULL;
	cudaCheckError(cudaMalloc((void**)&wave->ubuf_d, 3*(wave->nx)*(wave->ny)*sizeof(Number_t)));
	wave->isBulk_d = NULL;
	cudaCheckError(cudaMalloc((void**)&wave->isBulk_d, (wave->nx)*(wave->ny)*sizeof(bool)));
	wave->gradient_d = NULL;
	cudaCheckError(cudaMalloc((void**)&wave->gradient_d, 2*(wave->nx)*(wave->ny)*sizeof(Number_t)));

	//Set the pressures
	for(int i = 0; i < ny; i++){
		Number_t y = wave_sim_get_y(wave, i);
		for(int j = 0; j < nx; j++){
			Number_t x = wave_sim_get_x(wave, j);
			double dist = (x-0.5)*(x-0.5)+(y-0.5)*(y-0.5);
			if(dist < radius*radius && dist > radius*radius/2 && y < 0.7){
				wave->isBulk[j + nx*i] = false;
			} else{
				wave->isBulk[j + nx*i] = true;
			}
			wave->gradient[2*(j + nx*i)] = gradient(x+wave->dx/2, y, 0);
			wave->gradient[2*(j + nx*i)+1] = gradient(x, y+wave->dy/2, 1);
		}
	}

	Number_t * u = wave->ubuf;
	for(int i = 0; i < ny; i++){
		Number_t y = wave_sim_get_y(wave, i);
		for(int j = 0; j < nx; j++){
			Number_t x = wave_sim_get_x(wave, j);
			u[3*(j + i*nx)] = init_function(x, y, ctx);
		}
	}
	//Set the velocities
	for(int i = 0; i < ny; i++){
		for(int j = 0; j < nx; j++){
			u[3*(j + i*nx) + 1] = 0;
			u[3*(j + i*nx) + 2] = 0;
		}
	}

	wave_sim_apply_boundary(wave);

	Number_t * u_d = wave->ubuf_d;
	cudaCheckError(cudaMemcpy(u_d, u, 3*(wave->nx)*(wave->ny)*sizeof(Number_t), cudaMemcpyHostToDevice));
	cudaCheckError(cudaMemcpy(wave->isBulk_d, wave->isBulk, (wave->nx)*(wave->ny)*sizeof(bool), cudaMemcpyHostToDevice));
	cudaCheckError(cudaMemcpy(wave->gradient_d, wave->gradient, 2*(wave->nx)*(wave->ny)*sizeof(Number_t), cudaMemcpyHostToDevice));

	Number_t cco[11] = {wave->c, wave->dt,
				 	wave->dx, wave->dy,
				 	wave->xmin, wave->xmax,
				 	wave->ymin, wave->ymax,
				 	wave->pml_strength,
				 	wave->pml_width,
				 	wave->density};

	cudaMemcpyToSymbol(kernel_constants, cco, 11*sizeof(Number_t));

	return wave;
}

void wave_sim_free(Cuda_PML_Wave_2d_t wave){
	cudaFreeHost(wave->ubuf);
	cudaFree(wave->ubuf_d);
	free(wave);
}

Number_t wave_sim_get_x(Cuda_PML_Wave_2d_t wave, int j){
	return ((j*wave->xmax + (wave->nx - j)*wave->xmin)/wave->nx) + wave->dx/2;
}

Number_t wave_sim_get_y(Cuda_PML_Wave_2d_t wave, int i){
	return ((i*wave->ymax + (wave->ny - i)*wave->ymin)/wave->ny) + wave->dy/2;
}

void wave_sim_step(Cuda_PML_Wave_2d_t wave){
	//Copy to GPU
	Number_t * u = wave->ubuf;
	Number_t * u_d = wave->ubuf_d;

	cudaCheckError(cudaMemcpy(u_d, u, 3*(wave->nx)*(wave->ny)*sizeof(Number_t), cudaMemcpyHostToDevice));
	
	size_t blocks_x = ceil(wave->nx/16.0);
	size_t blocks_y = ceil(wave->ny/16.0);
	dim3 gridDim(blocks_x, blocks_y, 1);
	size_t threads_x = 16;
	size_t threads_y = 16;

	dim3 blockDim(threads_x, threads_y, 1);
	cuda_pml_wave_2d_velocity_kernel<<< gridDim, blockDim >>>(wave->ubuf_d,
															  wave->gradient_d,
															  wave->isBulk_d,
															  wave->t,
													 		  wave->nx, wave->ny);
	cudaCheckError(cudaGetLastError());
	cuda_pml_wave_2d_pressure_kernel<<< gridDim, blockDim >>>(wave->ubuf_d,
															  wave->isBulk_d,
													 wave->nx, wave->ny);
	cudaCheckError(cudaGetLastError());

	//Copy back
	cudaCheckError(cudaMemcpy(u, u_d, 3*(wave->nx)*(wave->ny)*sizeof(Number_t), cudaMemcpyDeviceToHost));
	
	wave->t += wave->dt;
	wave_sim_apply_boundary(wave);
	
}

void wave_sim_apply_boundary(Cuda_PML_Wave_2d_t wave){
	Number_t * u = wave->ubuf;
	int nx = wave->nx;
	int ny = wave->ny;

	for(int i = 0; i < ny; i++){
		u[3*(nx-1 + nx*i)+1] = 0;
	}

	for(int j = 0; j < nx; j++){
		u[3*(j + nx*(ny-1))+2] = 0;
	}
}