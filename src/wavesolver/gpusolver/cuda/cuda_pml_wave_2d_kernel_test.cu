#include "cuda_pml_wave_2d_kernel_math.cu"

__device__ __forceinline__ Number_t w_get_pos(int j, int nn, Number_t vmin, Number_t vmax, Number_t dd){
	return ((j*vmax + (nn - j)*vmin)/nn) + dd/2;
}

__constant__ Number_t kernel_constants[12];

__global__ void cuda_pml_wave_2d_velocity_kernel(Number_t * u,
												 Number_t * grad,
												 bool * isBulk,
												 Number_t t,
												 const int nx, const int ny){

	//__shared__ Number_t cache[BDIMX + 2][BDIMY + 2];
	Number_t local[2][2];

	Number_t dt = kernel_constants[1];
	Number_t idt = 1/dt;
	Number_t dx = kernel_constants[2];
	Number_t dy = kernel_constants[3];
	Number_t xmin = kernel_constants[4];
	Number_t xmax = kernel_constants[5];
	Number_t ymin = kernel_constants[6];
	Number_t ymax = kernel_constants[7];
	Number_t pml_strength = kernel_constants[8];
	Number_t pml_width = kernel_constants[9];
	Number_t density = kernel_constants[10];

	int j = blockIdx.x*blockDim.x + threadIdx.x;
	int i = blockIdx.y*blockDim.y + threadIdx.y;

	if(j < nx && i < ny){

		int idx = 3*(j + nx*i);

		//Set the position
		Number_t bx = w_get_pos(j, nx, xmin, xmax, dx);
		Number_t by = w_get_pos(i, ny, ymin, ymax, dy);

		local[0][0] = u[idx];
		Number_t oldVx = u[idx+1];
		Number_t oldVy = u[idx+2];

		bool ibulk = isBulk[j + i*nx];
		//Update velocities
		{
			Number_t newVx = 0;
			Number_t newVy = 0;

			Number_t absortion;
			Number_t update;
			Number_t gradient;
			//X
			if(j != nx-1){
				bool otherbulk = isBulk[(j+1)+i*nx];
				if(ibulk && otherbulk){
					local[1][0] = u[idx+3];
					absortion = pml_wave_2d_absortion(bx+dx/2, xmin, xmax, pml_strength, pml_width);
					update = pml_wave_2d_vel_update(idt, absortion);
					gradient = pml_wave_2d_gradient(idt, absortion, dx, density);

					newVx = oldVx*update + gradient*(local[1][0]-local[0][0]);
				} else if(ibulk){
					if(t < 3.14/100){
						newVx = -10*sin(100*t)*grad[2*(j + i*nx)]*dt + oldVx;
					}
				} else if(otherbulk){
					if(t < 3.14/100){
						newVx = 10*sin(100*t)*grad[2*(j + i*nx)]*dt + oldVx;
					}
				}
				u[idx+1] = newVx;
			}

			//Y
			if(i != ny-1){
				bool otherbulk = isBulk[j+(i+1)*nx];
				if(ibulk && otherbulk){
					local[0][1] = u[idx + 3*nx];
					absortion = pml_wave_2d_absortion(by+dy/2, ymin, ymax, pml_strength, pml_width);
					update = pml_wave_2d_vel_update(idt, absortion);
					gradient = pml_wave_2d_gradient(idt, absortion, dy, density);

					newVy = oldVy*update + gradient*(local[0][1]-local[0][0]);
				}else if(ibulk){
					newVy = 0;
				}else if(otherbulk){
					newVy = 0;
				}
				u[idx+2] = newVy;
			}
		}
	}
}

__global__ void cuda_pml_wave_2d_pressure_kernel(Number_t * u,
												 bool * isBulk,
												 const int nx, const int ny){

	//__shared__ Number_t cache[BDIMX + 2][BDIMY + 2];
	Number_t local[2][2];

	Number_t c = kernel_constants[0];
	Number_t dt = kernel_constants[1];
	Number_t idt = 1/dt;
	Number_t dx = kernel_constants[2];
	Number_t dy = kernel_constants[3];
	Number_t xmin = kernel_constants[4];
	Number_t xmax = kernel_constants[5];
	Number_t ymin = kernel_constants[6];
	Number_t ymax = kernel_constants[7];
	Number_t pml_strength = kernel_constants[8];
	Number_t pml_width = kernel_constants[9];
	Number_t density = kernel_constants[10];

	int j = blockIdx.x*blockDim.x + threadIdx.x;
	int i = blockIdx.y*blockDim.y + threadIdx.y;

	if(j < nx && i < ny && isBulk[j + i*nx]){

		int idx = 3*(j + nx*i);

		Number_t update = 0;
		Number_t temp = 0;

		Number_t oldU = u[idx];
		local[0][1] = u[idx+1];
		local[1][1] = u[idx+2];
		if(j != 0){
			local[0][0] = u[idx-3+1];
		} else{
			local[0][0] = 0;
		}

		//Set the position
		Number_t abs_d;
		Number_t dir_d;
		Number_t upd_d;
		Number_t div_d;
		//Update pressure
		{
			Number_t bx = w_get_pos(j, nx, xmin, xmax, dx);
			abs_d = pml_wave_2d_absortion(bx+dx/2, xmin, xmax, pml_strength, pml_width);
			dir_d = pml_wave_2d_directional(idt, abs_d);
			upd_d = pml_wave_2d_pre_update(idt, abs_d, dir_d);
			div_d = pml_wave_2d_pre_divergence(density, c, dir_d, dx);
			update += upd_d/2;
			temp += div_d*(local[0][1] - local[0][0]);
		}
		if(i != 0){
			local[1][0] = u[idx - 3*nx + 2];
		} else{
			local[1][0] = 0;
		}
		{
			Number_t by = w_get_pos(i, ny, ymin, ymax, dy);
			abs_d = pml_wave_2d_absortion(by+dy/2, ymin, ymax, pml_strength, pml_width);
			dir_d = pml_wave_2d_directional(idt, abs_d);
			upd_d = pml_wave_2d_pre_update(idt, abs_d, dir_d);
			div_d = pml_wave_2d_pre_divergence(density, c, dir_d, dy);
			update += upd_d/2;
			temp += div_d*(local[1][1] - local[1][0]);
		}

		//Write back to the global memory
		// u[idx] = threadIdx.x;
		Number_t newU = oldU*update + temp;
		u[ idx ] = newU;
	}
}