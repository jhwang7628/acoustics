#include "cuda_PAT_wave_3d.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cstdio>
#include <complex>
#include <boost/math/special_functions/bessel.hpp>

#include "cuda_helper.h"

#include "cuda_PAT_wave_3d_kernel.cu"

struct Cuda_PAT_Wave_3d_sim_data_t {
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

	Number_t * amplitude;
	Number_t * amplitude_d;

	Number_t * phase;
	Number_t * phase_d;

	bool * isBulk;
	bool * isBulk_d;

	Number_t * gradient;
	Number_t * gradient_d;

	int listening_count;
	Number_t * listening_positions_d;
	Number_t * listeningOutput;
	Number_t * listeningOutput_d;

	int num_multipole_coef;
	Number_t * multipole_coef;
	Number_t * multipole_coef_d;
	Number_t * integral_multipole_d;
	Number_t multipole_radius;


	Number_t xcenter, ycenter, zcenter;
	Number_t frequency;

	bool updated;
};

Number_t * wave_sim_get_u(Cuda_PAT_Wave_3d_t wave){
	if(wave->updated){
		return wave->ubuf;
	} else{
		
		cudaCheckError(cudaMemcpy(wave->ubuf, wave->ubuf_d, 4*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyDeviceToHost));
		cudaCheckError(cudaMemcpy(wave->phase, wave->phase_d, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyDeviceToHost));
		cudaCheckError(cudaMemcpy(wave->amplitude, wave->amplitude_d, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyDeviceToHost));
		wave->updated = true;

		return wave->ubuf;
	}
}

Number_t * wave_sim_get_amplitudes(Cuda_PAT_Wave_3d_t wave){
	if(wave->updated){
		return wave->amplitude;
	} else{
		
		cudaCheckError(cudaMemcpy(wave->ubuf, wave->ubuf_d, 4*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyDeviceToHost));
		cudaCheckError(cudaMemcpy(wave->phase, wave->phase_d, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyDeviceToHost));
		cudaCheckError(cudaMemcpy(wave->amplitude, wave->amplitude_d, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyDeviceToHost));
		wave->updated = true;

		return wave->amplitude;
	}	
}

Number_t * wave_sim_get_phases(Cuda_PAT_Wave_3d_t wave){
	if(wave->updated){
		return wave->phase;
	} else{
		
		cudaCheckError(cudaMemcpy(wave->ubuf, wave->ubuf_d, 4*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyDeviceToHost));
		cudaCheckError(cudaMemcpy(wave->phase, wave->phase_d, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyDeviceToHost));
		cudaCheckError(cudaMemcpy(wave->amplitude, wave->amplitude_d, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyDeviceToHost));
		wave->updated = true;

		return wave->phase;
	}	
}

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
								 Number_t multipole_radius){
	Cuda_PAT_Wave_3d_t wave = (Cuda_PAT_Wave_3d_t) malloc(sizeof(Cuda_PAT_Wave_3d_sim_data_t));
	
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
	wave->density = 1.22521;

	wave->nx = ceil((xmax-xmin)/cellsize);
	wave->ny = ceil((ymax-ymin)/cellsize);
	wave->nz = ceil((zmax-zmin)/cellsize);
	int nx = wave->nx;
	int ny = wave->ny;
	int nz = wave->nz;
	printf(">> %d %d %d %lf\n", nx, ny, nz, cellsize);

	wave->dx = (xmax-xmin)/wave->nx;
	wave->dy = (ymax-ymin)/wave->ny;
	wave->dz = (zmax-zmin)/wave->nz;

	wave->pml_width = pml_width;
	wave->pml_strength = pml_strength;
	wave->frequency = frequency;

	wave->ubuf = NULL;
	cudaCheckError(cudaMallocHost((void**)&wave->ubuf, 4*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t)));
	wave->ubuf_d = NULL;
	cudaCheckError(cudaMalloc((void**)&wave->ubuf_d, 4*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t)));

	wave->amplitude = NULL;
	cudaCheckError(cudaMallocHost((void**)&wave->amplitude, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t)));
	wave->amplitude_d = NULL;
	cudaCheckError(cudaMalloc((void**)&wave->amplitude_d, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t)));

	wave->phase = NULL;
	cudaCheckError(cudaMallocHost((void**)&wave->phase, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t)));
	wave->phase_d = NULL;
	cudaCheckError(cudaMalloc((void**)&wave->phase_d, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t)));

	wave->isBulk = NULL;
	cudaCheckError(cudaMallocHost((void**)&wave->isBulk, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(bool)));
	wave->isBulk_d = NULL;
	cudaCheckError(cudaMalloc((void**)&wave->isBulk_d, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(bool)));

	wave->gradient = NULL;
	cudaCheckError(cudaMallocHost((void**)&wave->gradient, 3*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t)));
	wave->gradient_d = NULL;
	cudaCheckError(cudaMalloc((void**)&wave->gradient_d, 3*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t)));

	wave->listening_count = listening_count;
	wave->listening_positions_d = NULL;
	wave->listeningOutput = NULL;
	wave->listeningOutput_d = NULL;
	if(listening_count > 0){
		cudaCheckError(cudaMalloc((void**)&wave->listening_positions_d, 3*listening_count*sizeof(Number_t)));
		cudaCheckError(cudaMallocHost((void**)&wave->listeningOutput, listening_count*sizeof(Number_t)));
		cudaCheckError(cudaMalloc((void**)&wave->listeningOutput_d, listening_count*sizeof(Number_t)));
	}

	wave->num_multipole_coef = num_multipole_coef;
	wave->multipole_coef = NULL;
	wave->multipole_coef_d = NULL;
	wave->integral_multipole_d = NULL;
	wave->multipole_radius = multipole_radius;
	if(num_multipole_coef > 0){
		int nn = (num_multipole_coef+1)*(num_multipole_coef+1);
		cudaCheckError(cudaMallocHost((void**)&wave->multipole_coef, 2*nn*sizeof(Number_t)));
		cudaCheckError(cudaMalloc((void**)&wave->multipole_coef_d, 2*nn*sizeof(Number_t)));
		cudaCheckError(cudaMalloc((void**)&wave->integral_multipole_d, 2*nn*wave->nx*sizeof(Number_t)));
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

	memset(u, 0, 4*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t));
	memset(wave->amplitude, 0, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t));
	memset(wave->phase, 0, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t));

	for(int k = 0; k < nz; k++){
		Number_t z = wave_sim_get_z(wave, k);
		for(int j = 0; j < ny; j++){
			Number_t y = wave_sim_get_y(wave, j);
			for(int i = 0; i < nx; i++){
				Number_t x = wave_sim_get_x(wave, i);
				Number_t val = initial(x, y, z);
				int idx = 4*(i + nx*(j + ny*k));
				u[idx] = val;
			}
		}
	}

	cudaCheckError(cudaMemcpy(wave->ubuf_d, wave->ubuf, 4*(wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyHostToDevice));
	cudaCheckError(cudaMemcpy(wave->amplitude_d, wave->amplitude, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyHostToDevice));
	cudaCheckError(cudaMemcpy(wave->phase_d, wave->phase, (wave->nx)*(wave->ny)*(wave->nz)*sizeof(Number_t), cudaMemcpyHostToDevice));
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
					 	wave->frequency,
					 	wave->dz};

	cudaMemcpyToSymbol(kernel_constants, cco, 18*sizeof(Number_t));

	return wave;
}

void wave_sim_free(Cuda_PAT_Wave_3d_t wave){
	cudaFreeHost(wave->ubuf);
	cudaFreeHost(wave->isBulk);
	cudaFreeHost(wave->gradient);
	cudaFreeHost(wave->listeningOutput);
	cudaFreeHost(wave->amplitude);
	cudaFreeHost(wave->phase);
	cudaFree(wave->ubuf_d);
	cudaFree(wave->isBulk_d);
	cudaFree(wave->gradient_d);
	cudaFree(wave->listeningOutput_d);
	cudaFree(wave->listening_positions_d);
	cudaFree(wave->amplitude_d);
	cudaFree(wave->phase_d);
		cudaFreeHost(wave->multipole_coef);
		cudaFree(wave->multipole_coef_d);
		cudaFree(wave->integral_multipole_d);
	free(wave);
}

Number_t wave_sim_get_x(Cuda_PAT_Wave_3d_t wave, int i){
	return ((i*wave->xmax + (wave->nx - i)*wave->xmin)/wave->nx) + wave->dx/2;
}

Number_t wave_sim_get_y(Cuda_PAT_Wave_3d_t wave, int j){
	return ((j*wave->ymax + (wave->ny - j)*wave->ymin)/wave->ny) + wave->dy/2;
}

Number_t wave_sim_get_z(Cuda_PAT_Wave_3d_t wave, int k){
	return ((k*wave->zmax + (wave->nz - k)*wave->zmin)/wave->nz) + wave->dz/2;
}


void wave_sim_step(Cuda_PAT_Wave_3d_t wave){
	size_t blocks_x = ceil(wave->nx/16.0);
	size_t blocks_y = ceil(wave->ny/16.0);
	dim3 gridDim(blocks_x, blocks_y, 1);
	size_t threads_x = 16;
	size_t threads_y = 16;

	printf("%lf -> %lf\n", sin(wave->t*wave->frequency), (wave->t/wave->dt) - (wave->multipole_radius/(wave->dt*wave->c)));
	dim3 blockDim(threads_x, threads_y, 1);
	cuda_pat_wave_3d_velocity_kernel<<< gridDim, blockDim >>>(wave->ubuf_d,
															  wave->gradient_d,
															  wave->isBulk_d,
															  wave->t,
													 		  wave->nx,
													 		  wave->ny,
													 		  wave->nz);
	cudaCheckError(cudaGetLastError());
	cuda_pat_wave_3d_pressure_kernel<<< gridDim, blockDim >>>(wave->ubuf_d,
															  wave->amplitude_d,
															  wave->phase_d,
															  wave->isBulk_d,
															  wave->nx,
															  wave->ny,
															  wave->nz,
															  wave->t,
															  (wave->t/wave->dt) - (wave->multipole_radius/(wave->dt*wave->c)));
	cudaCheckError(cudaGetLastError());

	wave->updated = false;
	wave->t += wave->dt;
}

void wave_set_listening_positions(Cuda_PAT_Wave_3d_t wave, int listening_count, Number_t * listening_positions){
	if(wave->listening_positions_d != NULL){
		cudaFree(wave->listening_positions_d);
		wave->listening_positions_d = NULL;
	}
	if(wave->listeningOutput_d != NULL){
		cudaFree(wave->listeningOutput_d);
		wave->listeningOutput_d = NULL;
	}
	if(wave->listeningOutput != NULL){
		cudaFreeHost(wave->listeningOutput);
		wave->listeningOutput = NULL;
	}

	wave->listening_count = listening_count;
	if(listening_count > 0){
		cudaCheckError(cudaMalloc((void**)&wave->listening_positions_d, 3*listening_count*sizeof(Number_t)));
		cudaCheckError(cudaMallocHost((void**)&wave->listeningOutput, listening_count*sizeof(Number_t)));
		cudaCheckError(cudaMalloc((void**)&wave->listeningOutput_d, listening_count*sizeof(Number_t)));
		cudaCheckError(cudaMemcpy(wave->listening_positions_d, listening_positions, 3*listening_count*sizeof(Number_t), cudaMemcpyHostToDevice));
	}
}

Number_t * wave_listen(Cuda_PAT_Wave_3d_t wave){
	if(wave->listening_count > 0){
		size_t blocks_x = ceil(wave->listening_count/256.0);
		dim3 gridDim(blocks_x, 1, 1);
		size_t threads_x = 256;

		dim3 blockDim(threads_x, 1, 1);

		cuda_pat_wave_3d_listen_kernel<<<gridDim, blockDim>>>(wave->ubuf_d,
															  wave->listeningOutput_d,
															  wave->listening_count,
															  wave->listening_positions_d,
															  wave->nx,
															  wave->ny,
															  wave->nz);
		cudaCheckError(cudaGetLastError());

		//Copy back
		cudaCheckError(cudaMemcpy(wave->listeningOutput, wave->listeningOutput_d, wave->listening_count*sizeof(Number_t), cudaMemcpyDeviceToHost));
		return wave->listeningOutput;
	} else{
		return NULL;
	}
}

void wave_sim_get_divisions(const Cuda_PAT_Wave_3d_t wave, int * nx, int * ny, int * nz){
	(*nx) = wave->nx;
	(*ny) = wave->ny;
	(*nz) = wave->nz;
}

Number_t wave_sim_get_current_time(const Cuda_PAT_Wave_3d_t wave){
	return wave->t;
}

#include <multipole/MultipoleMath.h>

using Multipole::spherical_harmonics;
using Multipole::hankel_1st;
using Multipole::hankel_2nd;
using Multipole::spherical_bessel;
using Multipole::regular_basis;
using Multipole::regular_basis_dir_deriv;

Number_t * wave_compute_multipole(Cuda_PAT_Wave_3d_t wave, Number_t radius){
	radius = wave->multipole_radius;
	//Compute the spherical bessel function here
	Number_t r = radius;
	Number_t ko = wave->frequency/wave->c;

	Number_t bessel_h[200];

	for(int l = 0; l < 160; l++){
		bessel_h[l] = spherical_bessel(l, r*ko);
	}

	cudaMemcpyToSymbol(bessel, bessel_h, 160*sizeof(Number_t));

	size_t blocks_x = ceil(wave->nx);
	size_t blocks_y = 1;
	dim3 gridDim(blocks_x, blocks_y, 1);
	size_t threads_x = wave->num_multipole_coef+1;
	size_t threads_y = 1;
	dim3 blockDim(threads_x, threads_y, 1);

	cuda_pat_wave_3d_compute_multipoles<<<gridDim, blockDim>>>(wave->integral_multipole_d,
															   wave->amplitude_d,
															   wave->phase_d,
															   wave->num_multipole_coef,
															   wave->nx,
															   wave->ny,
															   wave->nz,
															   radius);
	cudaCheckError(cudaGetLastError());

	int nn = (wave->num_multipole_coef + 1)*(wave->num_multipole_coef + 1);

	size_t blocks_x_2 = ceil(nn/32.0);
	size_t blocks_y_2 = 1;
	dim3 gridDim_2(blocks_x_2, blocks_y_2, 1);
	size_t threads_x_2 = 32;
	size_t threads_y_2 = 1;
	dim3 blockDim_2(threads_x_2, threads_y_2, 1);

	cuda_pat_wave_3d_reduce_multipoles<<<gridDim_2, blockDim_2>>>(wave->multipole_coef_d,
																  wave->integral_multipole_d,
																  wave->num_multipole_coef,
																  wave->nx);
	cudaCheckError(cudaGetLastError());

	cudaCheckError(cudaMemcpy(wave->multipole_coef, wave->multipole_coef_d, 2*nn*sizeof(Number_t), cudaMemcpyDeviceToHost));

	return wave->multipole_coef;
}

void wave_compute_contrib_slow(Cuda_PAT_Wave_3d_t wave, Number_t radius, Number_t x, Number_t y, Number_t z, int multipole_to_use){
	double cx = x - wave->xcenter;
	double cy = y - wave->ycenter;
	double cz = z - wave->zcenter;

	double r = sqrt(cx*cx + cy*cy + cz*cz);
	double theta = acos(cz/r);
	double phi = atan2(cy, cx);

	double ko = wave->frequency/wave->c;

	const int i = (int) floor((x - wave->dx/2 - wave->xmin)/wave->dx);
	const int j = (int) floor((y - wave->dy/2 - wave->ymin)/wave->dy);
	const int k = (int) floor((z - wave->dz/2 - wave->zmin)/wave->dz);

	double amp;
	double phase;
	std::complex<double> p;
	std::complex<double> dp_dn(0, 0);
	std::complex<double> temp;

	Number_t * amplitudes = wave_sim_get_amplitudes(wave);
	Number_t * phases = wave_sim_get_phases(wave);
	{
		amp = amplitudes[i + wave->nx*(j + wave->ny*k)];
		phase = phases[i + wave->nx*(j + wave->ny*k)];
		p = std::complex<double>(amp*cos(phase), amp*sin(phase));
	}

	{
		amp = amplitudes[i+1 + wave->nx*(j + wave->ny*k)];
		phase = phases[i+1 + wave->nx*(j + wave->ny*k)];
		temp = std::complex<double>(amp*cos(phase), amp*sin(phase));
		dp_dn += ((temp - p)*cx)/(r*wave->dx);
	}

	{
		amp = amplitudes[i + wave->nx*(j+1 + wave->ny*k)];
		phase = phases[i + wave->nx*(j+1 + wave->ny*k)];
		temp = std::complex<double>(amp*cos(phase), amp*sin(phase));
		dp_dn += ((temp - p)*cy)/(r*wave->dy);
	}

	{
		amp = amplitudes[i + wave->nx*(j + wave->ny*(k+1))];
		phase = phases[i + wave->nx*(j + wave->ny*(k+1))];
		temp = std::complex<double>(amp*cos(phase), amp*sin(phase));
		dp_dn += ((temp - p)*cz)/(r*wave->dz);
	}
	const double area = (abs(wave->dx * wave->dy * cz) + abs(wave->dx * cy *wave->dz) + abs(cx * wave->dy *wave->dz))/r;
	for(int l = 0; l <= multipole_to_use; l++){
		for(int m = -l; m <= l; m++){
			std::complex<double> contrib
			= ko*area*(
				regular_basis(-m, l, ko, r, theta, phi) * dp_dn
				- regular_basis_dir_deriv(-m, l, ko, r, theta, phi, cx/r, cy/r, cz/r)*p);
			// int mm = -m;
			
			// std::complex<double> contrib
			// = ko*(
			// 	regular_basis(-m, l, ko, r, theta, phi));
			

			// std::complex<double> contrib
			// = -ko*(regular_basis_dir_deriv(-m, l, ko, r, theta, phi, cx/r, cy/r, cz/r));
		

			wave->multipole_coef[2*(l*(l+1) + m)] += -contrib.imag();
			wave->multipole_coef[2*(l*(l+1) + m) + 1] += contrib.real();
		}
	}	
}

Number_t * wave_compute_multipole_slow(Cuda_PAT_Wave_3d_t wave, Number_t radius, int multipole_to_use){
	for(int i = 0; i < 2*(multipole_to_use+1)*(multipole_to_use+1); i++){
		wave->multipole_coef[i] = 0;
	}
	for(int k = 0; k < wave->nz; k++){
		for(int i = 0; i < wave->nx; i++){
			Number_t bx = wave_sim_get_x(wave, i);
			Number_t bz = wave_sim_get_z(wave, k);

			Number_t cx = bx - wave->xcenter;
			Number_t cz = bz - wave->zcenter;

			const Number_t sqdist = radius*radius - (cx)*(cx) - (cz)*(cz);
			if(sqdist > 0){
				const Number_t ya = wave->ycenter + sqrt(sqdist);
				wave_compute_contrib_slow(wave, radius, bx, ya, bz, multipole_to_use);
				const Number_t yb = wave->ycenter - sqrt(sqdist);
				wave_compute_contrib_slow(wave, radius, bx, yb, bz, multipole_to_use);
			} else if(sqdist == 0){
				const Number_t ya = wave->ycenter;
				wave_compute_contrib_slow(wave, radius, bx, ya, bz, multipole_to_use);
			}
		}
	}
	return wave->multipole_coef;
}

void wave_estimate_ijk(const Cuda_PAT_Wave_3d_t wave, Number_t x, Number_t y, Number_t z, int * i, int * j, int * k){
	*i = max(0, min((int) floor((x - wave->dx/2 - wave->xmin)/wave->dx), wave->nx));
	*j = max(0, min((int) floor((y - wave->dy/2 - wave->ymin)/wave->dy), wave->ny));
	*k = max(0, min((int) floor((z - wave->dz/2 - wave->zmin)/wave->dz), wave->nz));
}

void wave_estimate_with_multipole(Cuda_PAT_Wave_3d_t wave, Number_t x, Number_t y, Number_t z, Number_t * amplitude, Number_t * phase){
	double radius = wave->multipole_radius;
	Number_t * coef = wave_compute_multipole(wave, radius);
	int multipole_to_use = wave->num_multipole_coef;
	double cx = x - wave->xcenter;
	double cy = y - wave->ycenter;
	double cz = z - wave->zcenter;

	double ko = wave->frequency/wave->c;

	double r = sqrt(cx*cx + cy*cy + cz*cz);
	double theta = acos(cz / r);
	double phi  = atan2(cy, cx);

	std::complex<double> res(0, 0);
	for(int l = 0; l <= multipole_to_use; l++){
		std::complex<double> hank = hankel_2nd(l, r*ko);
		for(int m = -l; m <= l; m++){
			std::complex<double> spher = spherical_harmonics(m, l, theta, phi);
			std::complex<double> spherical(spher.real(), spher.imag());
			double re = coef[2*(l*(l+1) + m)];
			double im = coef[2*(l*(l+1) + m)+1];
			std::complex<double> multipole(re, im);
			res += hank*spherical*multipole;
		}
	}

	*amplitude = abs(res);
	*phase = arg(res);
}


void wave_test_multipole(Cuda_PAT_Wave_3d_t wave){
	Number_t radius = wave->multipole_radius;
	Number_t ko = wave->frequency/wave->c;
	int multipole_to_use = wave->num_multipole_coef;
	printf("Computing coefficients with radius %f and k %f\n", radius, ko);
	
	Number_t * coef;
	static bool use_gpu = true;
	if(use_gpu){
		coef = wave_compute_multipole(wave, radius);
	} else{
		coef = wave_compute_multipole_slow(wave, radius, multipole_to_use);
	}
	use_gpu = !use_gpu;
	printf("Evaluating multipole_coef\n");
	//Compute the spherical bessel function here
	for(int l = 0; l <= multipole_to_use; l++){
    				//printf("r = %f and k = %f\n", r, ko);
    				//printf("std::complex<double> hhank = boost::math::sph_hankel_2(%f, %f);\n", l, r*ko);
		for(int m = -l; m <= l; m++){
			Number_t re = coef[2*(l*(l+1) + m)];
			Number_t im = coef[2*(l*(l+1) + m)+1];
			printf("M(%d,%d) = %f + i%f\n", l, m, re, im);
		}
	}

	printf("USING %d\n", multipole_to_use);

	for(int k = 0; k < wave->nz; k++){
		if(k%50 == 0)printf("k: %d\n", k);
		for(int j = wave->ny/2-1; j < wave->ny/2+2; j++){
			for(int i = 0; i < wave->nx; i++){
				double bx = wave_sim_get_x(wave, i);
				double by = wave_sim_get_y(wave, j);
				double bz = wave_sim_get_z(wave, k);

				double cx = bx - wave->xcenter;
				double cy = by - wave->ycenter;
				double cz = bz - wave->zcenter;

				double r = sqrt(cx*cx + cy*cy + cz*cz);
				double theta = acos(cz / r);
    			double phi  = atan2(cy, cx);
    			if(r > wave->multipole_radius){
	    			std::complex<double> res(0, 0);
	    			for(int l = 0; l <= multipole_to_use; l++){
	    				//std::complex<double> hank(boost::math::sph_bessel(l, r*ko), boost::math::sph_neumann(l, r*ko));
	    				std::complex<double> hank = hankel_2nd(l, r*ko);
	    				for(int m = -l; m <= l; m++){
	    					std::complex<double> spher = spherical_harmonics(m, l, theta, phi);
	    					std::complex<double> spherical(spher.real(), spher.imag());
	    					double re = coef[2*(l*(l+1) + m)];
	    					double im = coef[2*(l*(l+1) + m)+1];
	    					std::complex<double> multipole(re, im);
	    					res += hank*spherical*multipole;
	    				}
	    			}
	    			wave->ubuf[4*(i + wave->nx*(j + wave->ny*k))] = res.real();
	    		} else if(r > 0.1*wave->multipole_radius){
	    			std::complex<double> res(0, 0);
	    			for(int l = 0; l <= multipole_to_use; l++){
	    				//std::complex<double> hank(boost::math::sph_bessel(l, r*ko), boost::math::sph_neumann(l, r*ko));
	    				std::complex<double> hank = hankel_2nd(l, r*ko);
	    				for(int m = -l; m <= l; m++){
	    					std::complex<double> spher = spherical_harmonics(m, l, theta, phi);
	    					std::complex<double> spherical(spher.real(), spher.imag());
	    					double re = coef[2*(l*(l+1) + m)];
	    					double im = coef[2*(l*(l+1) + m)+1];
	    					std::complex<double> multipole(re, im);
	    					res += hank*spherical*multipole;
	    				}
	    			}
	    			wave->ubuf[4*(i + wave->nx*(j + wave->ny*k))] = res.real();	    			
	    		}
			}
		}
	}
}

void wave_sim_get_bounds(const Cuda_PAT_Wave_3d_t wave,
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