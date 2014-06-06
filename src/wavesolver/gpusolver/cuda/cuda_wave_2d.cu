#include "cuda_wave_2d.h"

#include <cuda_runtime.h>

#include <ctime>
#include <cstdio>

#include "cuda_wave_2d_kernel.cu"

struct Cuda_Wave_2d_sim_data_t {
	Number_t xmin, ymin;
	Number_t xmax, ymax;

	Number_t dt;
	Number_t dx;
	Number_t dy;

	Number_t t;

	Number_t c;

	int nx, ny;

	Number_t * ubuf;

	int step;

	//CUDA
	Number_t * ubuf_d;

};


Number_t * wave_sim_get_u_device(Cuda_Wave_2d_t wave, int offset){
	return wave->ubuf_d + ((wave->step + offset + 3)%3)*(wave->nx+2)*(wave->ny+2);
}

Cuda_Wave_2d_t wave_sim_init(Number_t xmin, Number_t ymin,
						Number_t xmax, Number_t ymax,
						Number_t c, Number_t dt, 
						int nx, int ny,
						Number_t (*init_function)(Number_t, Number_t, void *),
						void * ctx){
	Cuda_Wave_2d_t wave = (Cuda_Wave_2d_t) malloc(sizeof(Cuda_Wave_2d_sim_data_t));
	
	wave->xmin = xmin;
	wave->ymin = ymin;
	wave->xmax = xmax;
	wave->ymax = ymax;

	wave->c = c;

	wave->dt = dt;
	wave->step = 0;

	wave->t = 0;

	wave->nx = nx;
	wave->ny = ny;

	wave->dx = (xmax-xmin)/nx;
	wave->dy = (ymax-ymin)/ny;
	int err;
	wave->ubuf = NULL;
	err = cudaMallocHost((void**)&wave->ubuf, 3*(wave->nx+2)*(wave->ny+2)*sizeof(Number_t));
	wave->ubuf_d = NULL;
	err = cudaMalloc((void**)&wave->ubuf_d, 3*(wave->nx+2)*(wave->ny+2)*sizeof(Number_t));


	Number_t * u = wave_sim_get_u(wave, 0);
	for(int i = 1; i <= ny; i++){
		Number_t y = wave_sim_get_y(wave, i);
		for(int j = 1; j <= nx; j++){
			Number_t x = wave_sim_get_x(wave, j);
			u[j + i*(nx+2)] = init_function(x, y, ctx);
		}
	}

	wave_sim_apply_boundary(wave);

	Number_t * uold_d = wave_sim_get_u_device(wave, -1);
	Number_t * u_d 	= wave_sim_get_u_device(wave,  0);
	Number_t * unew_d = wave_sim_get_u_device(wave,  1);

	cudaMemcpy(u_d, u, (wave->nx+2)*(wave->ny+2)*sizeof(Number_t), cudaMemcpyHostToDevice);
	cudaMemcpy(uold_d, u_d, (wave->nx+2)*(wave->ny+2)*sizeof(Number_t), cudaMemcpyDeviceToDevice);

	return wave;
}

void wave_sim_free(Cuda_Wave_2d_t wave){
	cudaFreeHost(wave->ubuf);
	cudaFree(wave->ubuf_d);
	free(wave);
}

Number_t wave_sim_get_x(Cuda_Wave_2d_t wave, int j){
	return ((j-1)*wave->xmax + (wave->nx - j)*wave->xmin)/(wave->nx-1);
}

Number_t wave_sim_get_y(Cuda_Wave_2d_t wave, int i){
	return ((i-1)*wave->ymax + (wave->ny - i)*wave->ymin)/(wave->ny-1);
}

Number_t * wave_sim_get_u(Cuda_Wave_2d_t wave, int offset){
	return wave->ubuf + ((wave->step + offset + 3)%3)*(wave->nx+2)*(wave->ny+2);
}



void wave_sim_step(Cuda_Wave_2d_t wave){
	Number_t * u = wave_sim_get_u(wave,  0);
	Number_t * unew = wave_sim_get_u(wave,  1);
	Number_t * uold_d = wave_sim_get_u_device(wave, -1);
	Number_t * u_d 	= wave_sim_get_u_device(wave,  0);
	Number_t * unew_d = wave_sim_get_u_device(wave,  1);
		//Copy to GPU
		int err = cudaMemcpy(u_d, u, (wave->nx+2)*(wave->ny+2)*sizeof(Number_t), cudaMemcpyHostToDevice);
		//Run kernel
		size_t blocks_x = ceil(wave->nx/16.0);
		size_t blocks_y = ceil(wave->ny/16.0);
		dim3 gridDim(blocks_x, blocks_y, 1);
		size_t threads_x = ceil(wave->nx/(float) blocks_x);
		size_t threads_y = ceil(wave->ny/(float) blocks_y);
		dim3 blockDim(threads_x, threads_y, 1);
		cuda_wave_2d_kernel<<< gridDim, blockDim >>>(unew_d, u_d, uold_d, wave->nx, wave->ny, wave->c, wave->dt, wave->dx, wave->dy);
		//Copy back
		err = cudaMemcpy(unew, unew_d, (wave->nx+2)*(wave->ny+2)*sizeof(Number_t), cudaMemcpyDeviceToHost);
	
		wave->step++;
		wave->t += wave->dt;
		wave_sim_apply_boundary(wave);
		
}

void wave_sim_apply_boundary(Cuda_Wave_2d_t wave){
	Number_t * u = wave_sim_get_u(wave,  0);
	int nx = wave->nx;
	int ny = wave->ny;

	for(int i = 1; i <= ny; i++){
		u[0 + i*(nx+2)] = u[2+i*(nx+2)];
		u[(nx+1) + i*(nx+2)] = u[(nx-1)+i*(nx+2)];
	}

	for(int j = 1; j <= nx; j++){
		u[j + 0*(nx+2)] = u[j+2*(nx+2)];
		u[j + (ny+1)*(nx+2)] = u[j+(ny-1)*(nx+2)];
	}
}