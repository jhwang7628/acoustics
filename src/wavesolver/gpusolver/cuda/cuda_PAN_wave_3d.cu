#include "cuda_PAN_wave_3d.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cstdio>

#include "cuda_helper.h"

#include "cuda_PAN_wave_3d_kernel.cu"

struct Cuda_PAN_Wave_3d_sim_data_t {
	Number_t xmin, ymin, zmin;
	Number_t xmax, ymax, zmax;

	Number_t dt;

	Number_t dx, dy, dz;

	Number_t t;

	Number_t c;

	int nx, ny, nz;

	Number_t pml_width;
	Number_t pml_strength;
	Number_t density;

	Number_t * ubuf;
	Number_t * ubuf_d;

	bool * isBulk;
	bool * isBulk_d;

	Number_t * gradient;
	Number_t * gradient_d;

	int listening_count;
	Number_t * listening_positions_d;
	Number_t * listeningOutput;
	Number_t * listeningOutput_d;

	Number_t xcenter, ycenter, zcenter;
	Number_t pulse;

	bool updated;
};

Number_t * wave_sim_get_u(Cuda_PAN_Wave_3d_t wave){
	if(wave->updated){
		return wave->ubuf;
	} else{
		Number_t * u = wave->ubuf;
		Number_t * u_d = wave->ubuf_d;

		cudaCheckError(cudaMemcpy(u, u_d, 4*6*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyDeviceToHost));
		wave->updated = true;
		return u;
	}
}

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
								 Number_t pulse){
	Cuda_PAN_Wave_3d_t wave = (Cuda_PAN_Wave_3d_t) malloc(sizeof(Cuda_PAN_Wave_3d_sim_data_t));
	
	wave->xmin = xmin;
	wave->ymin = ymin;
	wave->zmin = zmin;
	wave->xmax = xmax;
	wave->ymax = ymax;
	wave->zmax = zmax;

	wave->updated = true;

	wave->xcenter = xcenter;
	wave->ycenter = ycenter;
	wave->zcenter = zcenter;

	wave->c = c;

	wave->dt = dt;

	wave->t = 0;
	wave->density = 1;

	wave->nx = ceil((xmax-xmin)/cellsize);
	wave->ny = ceil((ymax-ymin)/cellsize);
	wave->nz = ceil((zmax-zmin)/cellsize);
	int nx = wave->nx;
	int ny = wave->ny;
	int nz = wave->nz;

	wave->dx = (xmax-xmin)/wave->nx;
	wave->dy = (ymax-ymin)/wave->ny;
	wave->dz = (zmax-zmin)/wave->nz;

	wave->pml_width = pml_width;
	wave->pml_strength = pml_strength;
	wave->pulse = pulse;

	wave->ubuf = NULL;
	cudaCheckError(cudaMallocHost((void**)&wave->ubuf, 4*6*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t)));
	wave->ubuf_d = NULL;
	cudaCheckError(cudaMalloc((void**)&wave->ubuf_d, 4*6*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t)));

	wave->isBulk = NULL;
	cudaCheckError(cudaMallocHost((void**)&wave->isBulk, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(bool)));
	wave->isBulk_d = NULL;
	cudaCheckError(cudaMalloc((void**)&wave->isBulk_d, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(bool)));

	wave->gradient = NULL;
	cudaCheckError(cudaMallocHost((void**)&wave->gradient, 3*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t)));
	wave->gradient_d = NULL;
	cudaCheckError(cudaMalloc((void**)&wave->gradient_d, 3*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t)));

	if(listening_count > 0){
		wave->listening_count = listening_count;

		wave->listening_positions_d = NULL;
		cudaCheckError(cudaMalloc((void**)&wave->listening_positions_d, 3*listening_count*sizeof(Number_t)));

		wave->listeningOutput = NULL;
		cudaCheckError(cudaMallocHost((void**)&wave->listeningOutput, 6*listening_count*sizeof(Number_t)));
		wave->listeningOutput_d = NULL;
		cudaCheckError(cudaMalloc((void**)&wave->listeningOutput_d, 6*listening_count*sizeof(Number_t)));
	}

	for(int k = 0; k < nz; k++){
		Number_t z = wave_sim_get_z(wave, k);
		for(int j = 0; j < ny; j++){
			Number_t y = wave_sim_get_y(wave, j);
			for(int i = 0; i < nx; i++){
				Number_t x = wave_sim_get_x(wave, i);
				if(boundary(x, y, z)){
					wave->isBulk[(i + nx*(j + ny*k))] = false;
				} else{
					wave->isBulk[(i + nx*(j + ny*k))] = true;
				}
				int idx = 3*(i + nx*(j + ny*k));
				wave->gradient[idx] = gradient(x+wave->dx/2, y, z, 0);
				wave->gradient[idx+1] = gradient(x, y+wave->dy/2, z, 1);
				wave->gradient[idx+2] = gradient(x, y, z+wave->dz/2, 2);
			}
		}
	}

	//Set the pressures
	Number_t * u = wave->ubuf;

	memset(u, 0, 4*6*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t));

	for(int k = 0; k < nz; k++){
		Number_t z = wave_sim_get_z(wave, k);
		for(int j = 0; j < ny; j++){
			Number_t y = wave_sim_get_y(wave, j);
			for(int i = 0; i < nx; i++){
				Number_t x = wave_sim_get_x(wave, i);
				Number_t val = initial(x, y, z);
				int idx = 4*(i + nx*(j + ny*k));
				int stride = 4*nx*ny*nz;
				u[idx] = val;
				u[idx+stride] = val;
				u[idx+2*stride] = val;
				u[idx+3*stride] = val;
				u[idx+4*stride] = val;
				u[idx+5*stride] = val;
			}
		}
	}

	cudaCheckError(cudaMemcpy(wave->ubuf_d, wave->ubuf, 4*6*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyHostToDevice));
	cudaCheckError(cudaMemcpy(wave->isBulk_d, wave->isBulk, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(bool), cudaMemcpyHostToDevice));
	cudaCheckError(cudaMemcpy(wave->gradient_d, wave->gradient, 3*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyHostToDevice));
	if(listening_count > 0){
		cudaCheckError(cudaMemcpy(wave->listening_positions_d, listening_positions, 3*listening_count*sizeof(Number_t), cudaMemcpyHostToDevice));	
	}
	
	Number_t cco[18] = {wave->c, wave->dt,
					 	wave->dx, wave->dy,
					 	wave->xmin, wave->xmax,
					 	wave->ymin, wave->ymax,
					 	wave->zmin, wave->zmax,
					 	wave->pml_strength,
					 	wave->pml_width,
					 	wave->density,
					 	wave->xcenter,
					 	wave->ycenter,
					 	wave->zcenter,
					 	wave->pulse,
					 	wave->dz};

	cudaMemcpyToSymbol(kernel_constants, cco, 18*sizeof(Number_t));

	return wave;
}

void wave_sim_free(Cuda_PAN_Wave_3d_t wave){
	cudaFreeHost(wave->ubuf);
	cudaFreeHost(wave->isBulk);
	cudaFreeHost(wave->gradient);
	cudaFreeHost(wave->listeningOutput);
	cudaFree(wave->ubuf_d);
	cudaFree(wave->isBulk_d);
	cudaFree(wave->gradient_d);
	cudaFree(wave->listeningOutput_d);
	cudaFree(wave->listening_positions_d);
	free(wave);
}

Number_t wave_sim_get_x(Cuda_PAN_Wave_3d_t wave, int i){
	return ((i*wave->xmax + (wave->nx - i)*wave->xmin)/wave->nx) + wave->dx/2;
}

Number_t wave_sim_get_y(Cuda_PAN_Wave_3d_t wave, int j){
	return ((j*wave->ymax + (wave->ny - j)*wave->ymin)/wave->ny) + wave->dy/2;
}

Number_t wave_sim_get_z(Cuda_PAN_Wave_3d_t wave, int k){
	return ((k*wave->zmax + (wave->nz - k)*wave->zmin)/wave->nz) + wave->dz/2;
}


void wave_sim_step(Cuda_PAN_Wave_3d_t wave){
	size_t blocks_x = ceil(wave->nx/16.0);
	size_t blocks_y = ceil(wave->ny/16.0);
	dim3 gridDim(blocks_x, blocks_y, 1);
	size_t threads_x = 16;
	size_t threads_y = 16;

	dim3 blockDim(threads_x, threads_y, 1);
	int stride = 4*(wave->nx)*(wave->ny)*(wave->nz);
	for(int i = 0; i < 6; i++){
		cuda_pan_wave_3d_velocity_kernel<<< gridDim, blockDim >>>(wave->ubuf_d+i*stride,
																  wave->gradient_d,
																  wave->isBulk_d,
																  wave->t,
														 		  wave->nx,
														 		  wave->ny,
														 		  wave->nz,
														 		  i);
		cudaCheckError(cudaGetLastError());
		cuda_pan_wave_3d_pressure_kernel<<< gridDim, blockDim >>>(wave->ubuf_d+i*stride,
																  wave->isBulk_d,
																  wave->nx,
																  wave->ny,
																  wave->nz);
	}
	cudaCheckError(cudaGetLastError());
	wave->updated = false;
	wave->t += wave->dt;
}

Number_t * wave_listen(Cuda_PAN_Wave_3d_t wave, int field){
	if(wave->listening_count > 0){
		size_t blocks_x = ceil(wave->listening_count/256.0);
		dim3 gridDim(blocks_x, 1, 1);
		size_t threads_x = 256;

		dim3 blockDim(threads_x, 1, 1);
		int stride = 4*(wave->nx)*(wave->ny)*(wave->nz)*field;
		int linstride = wave->listening_count*field;

		cuda_pan_wave_3d_listen_kernel<<<gridDim, blockDim>>>(wave->ubuf_d+stride,
															  wave->listeningOutput_d,
															  wave->listening_count,
															  wave->listening_positions_d,
															  wave->nx,
															  wave->ny,
															  wave->nz);
		cudaCheckError(cudaGetLastError());

		//Copy back
		cudaCheckError(cudaMemcpy(wave->listeningOutput+linstride, wave->listeningOutput_d, wave->listening_count*sizeof(Number_t), cudaMemcpyDeviceToHost));
		return wave->listeningOutput+linstride;
	} else{
		return NULL;
	}
}

void wave_sim_get_divisions(const Cuda_PAN_Wave_3d_t wave, int * nx, int * ny, int * nz){
	(*nx) = wave->nx;
	(*ny) = wave->ny;
	(*nz) = wave->nz;
}

Number_t wave_sim_get_current_time(const Cuda_PAN_Wave_3d_t wave){
	return wave->t;
}

void wave_sim_get_bounds(const Cuda_PAN_Wave_3d_t wave,
						 Number_t * xmin, Number_t * xmax,
						 Number_t * ymin, Number_t * ymax,
						 Number_t * zmin, Number_t * zmax){
	(*xmin) = wave->xmin;
	(*xmax) = wave->xmax;
	(*ymin) = wave->ymin;
	(*ymax) = wave->ymax;
	(*zmin) = wave->zmin;
	(*zmax) = wave->zmax;
}