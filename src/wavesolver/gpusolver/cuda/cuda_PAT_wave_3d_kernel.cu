__constant__ Number_t kernel_constants[18];
__constant__ Number_t bessel[200];
 
#include "cuda_PAT_wave_3d_kernel_math.cu"

__device__ __forceinline__ Number_t w_get_pos(int j, int nn, Number_t vmin, Number_t vmax, Number_t dd){
	return ((j*vmax + (nn - j)*vmin)/nn) + dd/2;
}


#define RRADX 16
#define RRADY 16

__global__ void cuda_pat_wave_3d_velocity_kernel(Number_t * __restrict__ u,
												 const Number_t * __restrict__ grad,
												 const bool * __restrict__ isBulk,
												 Number_t t,
												 const int nx,
												 const int ny,
												 const int nz){
	const Number_t dt = kernel_constants[1];
	const Number_t idt = 1/dt;
	const Number_t dx = kernel_constants[2];
	const Number_t dy = kernel_constants[3];
	const Number_t xmin = kernel_constants[4];
	const Number_t xmax = kernel_constants[5];
	const Number_t ymin = kernel_constants[6];
	const Number_t ymax = kernel_constants[7];
	const Number_t zmin = kernel_constants[8];
	const Number_t zmax = kernel_constants[9];
	const Number_t pml_strength = kernel_constants[10];
	const Number_t pml_width = kernel_constants[11];
	const Number_t density = kernel_constants[12];
	const Number_t frequency = kernel_constants[16];
	const Number_t dz = kernel_constants[17];
	const Number_t sint = frequency*cos(frequency*t);


	Number_t local_z[4];
	Number_t local_old[4];

	bool ibulk;
	bool bulk_z;

	__shared__ Number_t cache[RRADX+2][RRADY+2];
	__shared__ bool cache_bulk[RRADX+2][RRADY+2];

	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	const int j = blockIdx.y*blockDim.y + threadIdx.y;

	const int bdx = min(RRADX, nx-i);
	const int bdy = min(RRADY, ny-j);

	if(i < nx && j < ny){
		const Number_t bx = w_get_pos(i, nx, xmin, xmax, dx);
		const Number_t by = w_get_pos(j, ny, ymin, ymax, dy);

		const int tx = threadIdx.x + 1;
		const int ty = threadIdx.y + 1;

		const Number_t absortion_x = pat_wave_3d_absortion(bx+dx/2, xmin, xmax, pml_strength, pml_width);
		const Number_t update_x = pat_wave_3d_vel_update(idt, absortion_x);
		const Number_t gradient_x = pat_wave_3d_gradient(idt, absortion_x, dx, density);
		
		const Number_t absortion_y = pat_wave_3d_absortion(by+dy/2, ymin, ymax, pml_strength, pml_width);
		const Number_t update_y = pat_wave_3d_vel_update(idt, absortion_y);
		const Number_t gradient_y = pat_wave_3d_gradient(idt, absortion_y, dy, density);
		

		//Compute the first guy:
		{
			const int bbbase = 4*(i + nx*j);

			local_z[0] = u[bbbase+0];
			local_z[1] = u[bbbase+1];
			local_z[2] = u[bbbase+2];
			local_z[3] = u[bbbase+3];
			bulk_z = isBulk[(i + nx*j)];
		}

		#pragma unroll 2
		for(int other = 0; other < nz-1; other++){
			const int k = other;
			const int base = i + nx*(j + ny*k);
			const int idx = 4*base;
			const Number_t bz = w_get_pos(k, nz, zmin, zmax, dz);

			Number_t local_new[4];

			local_old[0] = local_z[0];
			local_old[1] = local_z[1];
			local_old[2] = local_z[2];
			local_old[3] = local_z[3];
			local_new[0] = local_old[0];
			local_new[1] = local_new[2] = local_new[3] = 0;

			ibulk = bulk_z;

			cache[tx][ty] = local_old[0];
			cache_bulk[tx][ty] = ibulk;

			__syncthreads();

			if(threadIdx.x == 0){
				if(i+bdx < nx){
					const int base = 4*((i+bdx) + nx*(j + ny*k));
					cache[tx+bdx][ty] = u[base+0];
					cache_bulk[tx+bdx][ty] = isBulk[(i+bdx) + nx*(j + ny*k)];
				} else{
					cache[tx+bdx][ty] = 0;
					cache_bulk[tx+bdx][ty] = false;
				}
			}
			if(threadIdx.y == 0){
				if(j+bdy < ny){
					const int base = 4*(i + nx*((j+bdy) + ny*k));
					cache[tx][ty+bdy] = u[base+0];
					cache_bulk[tx][ty+bdy] = isBulk[(i + nx*((j+bdy) + ny*k))];
				} else{
					cache[tx][ty+bdy] = 0;
					cache_bulk[tx][ty+bdy] = false;
				}
			}
			__syncthreads();

			//Solve for X
			if(i != nx-1){
				const bool otherbulk = cache_bulk[tx+1][ty];
				if(ibulk && otherbulk){
					local_new[1] = local_old[1]*update_x + gradient_x*(cache[tx+1][ty] - cache[tx][ty]);
				} else if(ibulk || otherbulk){
					Number_t gradi = grad[3*base]*sint;
					// if(ibulk){
					// 	gradi = -gradi;
					// }
					local_new[1] = gradi;
				}
			}
			//Solve for Y
			if(j != ny-1){
				const bool otherbulk = cache_bulk[tx][ty+1];
				if(ibulk && otherbulk){
					local_new[2] = local_old[2]*update_y + gradient_y*(cache[tx][ty+1] - cache[tx][ty]);
				} else if(ibulk || otherbulk){
					Number_t gradi = grad[3*base+1]*sint;
					// if(ibulk){
					// 	gradi = -gradi;
					// }
					local_new[2] = gradi;
				}
			}


			//Solve for Z
			{
				const int bbase = (i + nx*(j + ny*(k+1)));
				const int bidx = 4*bbase;
				local_z[0] = u[bidx+0];
				local_z[1] = u[bidx+1];
				local_z[2] = u[bidx+2];
				local_z[3] = u[bidx+3];
				bulk_z = isBulk[bbase];
				const bool otherbulk = bulk_z;
				if(ibulk && otherbulk){
					const Number_t absortion = pat_wave_3d_absortion(bz+dz/2, zmin, zmax, pml_strength, pml_width);
					const Number_t update = pat_wave_3d_vel_update(idt, absortion);
					const Number_t gradient = pat_wave_3d_gradient(idt, absortion, dz, density);
					local_new[3] = local_old[3]*update + gradient*(local_z[0] - cache[tx][ty]);
				} else if(ibulk || otherbulk){
					Number_t gradi = grad[3*base+2]*sint;
					// if(ibulk){
					// 	gradi = -gradi;
					// }
					local_new[3] = gradi;
				}
			}
			//u[idx + 0] = local_new[0];
			u[idx + 1] = local_new[1];
			u[idx + 2] = local_new[2];
			u[idx + 3] = local_new[3];
		}
	}
}

__device__ Number_t my_asin(Number_t x){
	x = ((x > 1) ? 1
				 : ((x < -1) ? -1
				 	         :  x
				 	         )
				 );
	return asin(x);
}


__global__ void cuda_pat_wave_3d_pressure_kernel(Number_t * __restrict__ u,
												 Number_t * __restrict__ amplitude,
												 Number_t * __restrict__ phase,
												 bool * __restrict__ isBulk,
												 const int nx,
												 const int ny,
												 const int nz,
												 Number_t t,
												 Number_t steps){

	const Number_t c = kernel_constants[0];
	const Number_t dt = kernel_constants[1];
	const Number_t idt = 1/dt;
	const Number_t dx = kernel_constants[2];
	const Number_t dy = kernel_constants[3];
	const Number_t xmin = kernel_constants[4];
	const Number_t xmax = kernel_constants[5];
	const Number_t ymin = kernel_constants[6];
	const Number_t ymax = kernel_constants[7];
	const Number_t zmin = kernel_constants[8];
	const Number_t zmax = kernel_constants[9];
	const Number_t pml_strength = kernel_constants[10];
	const Number_t pml_width = kernel_constants[11];
	const Number_t density = kernel_constants[12];
	const Number_t dz = kernel_constants[17];
	const Number_t frequency = kernel_constants[16];

	Number_t local_z;

	__shared__ Number_t cache[RRADX+2][RRADY+2][4];

	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	const int j = blockIdx.y*blockDim.y + threadIdx.y;

	const int bdx = min(RRADX, nx-i);
	const int bdy = min(RRADY, ny-j);

	if(i < nx && j < ny){
		const Number_t bx = w_get_pos(i, nx, xmin, xmax, dx);
		const Number_t abs_x = pat_wave_3d_absortion(bx+dx/2, xmin, xmax, pml_strength, pml_width);
		const Number_t dir_x = pat_wave_3d_directional(idt, abs_x);
		const Number_t upd_x = pat_wave_3d_pre_update(idt, abs_x, dir_x);
		const Number_t div_x = pat_wave_3d_pre_divergence(density, c, dir_x, dx);

		const Number_t by = w_get_pos(j, ny, ymin, ymax, dy);
		const Number_t abs_y = pat_wave_3d_absortion(by+dy/2, ymin, ymax, pml_strength, pml_width);
		const Number_t dir_y = pat_wave_3d_directional(idt, abs_y);
		const Number_t upd_y = pat_wave_3d_pre_update(idt, abs_y, dir_y);
		const Number_t div_y = pat_wave_3d_pre_divergence(density, c, dir_y, dy);

		const int tx = threadIdx.x + 1;
		const int ty = threadIdx.y + 1;

		{
			const int bbase = 4*(i + nx*j);
			cache[tx][ty][0] = u[bbase+0];
			cache[tx][ty][1] = u[bbase+1];
			cache[tx][ty][2] = u[bbase+2];
			cache[tx][ty][3] = u[bbase+3];
		}

		for(int other = 1; other < nz; other++){
			const int k = other;
			const int base = i + nx*(j+ny*k);
			const int idx = 4*base;

			{
				local_z = cache[tx][ty][3];
				cache[tx][ty][0] = u[idx+0];
				cache[tx][ty][1] = u[idx+1];
				cache[tx][ty][2] = u[idx+2];
				cache[tx][ty][3] = u[idx+3];
			}
			if(threadIdx.x == 0){
				if(i != 0){
					const int base = 4*((i-1) + nx*(j + ny*k));
					cache[0][ty][0] = u[base+0];
					cache[0][ty][1] = u[base+1];
					cache[0][ty][2] = u[base+2];
					cache[0][ty][3] = u[base+3];
				} else{
					cache[0][ty][0] = 0;
					cache[0][ty][1] = 0;
					cache[0][ty][2] = 0;
					cache[0][ty][3] = 0;
				}
			}
			if(threadIdx.y == 0){
				if(j != 0){
					const int base = 4*(i + nx*((j-1) + ny*k));
					cache[tx][0][0] = u[base+0];
					cache[tx][0][1] = u[base+1];
					cache[tx][0][2] = u[base+2];
					cache[tx][0][3] = u[base+3];
				} else{
					cache[tx][0][0] = 0;
					cache[tx][0][1] = 0;
					cache[tx][0][2] = 0;
					cache[tx][0][3] = 0;
				}
			}

			__syncthreads();
			const bool ibulk = isBulk[base];

			if(ibulk){
				Number_t update = 0;
				Number_t diver = 0;

				//Solve for X
				{
					diver += div_x*(cache[tx][ty][1]-cache[tx-1][ty][1]);
				}

				//Solve for Y
				{
					diver += div_y*(cache[tx][ty][2]-cache[tx][ty-1][2]);
				}

				//Solve for Z
				{
					const Number_t bz = w_get_pos(k, nz, zmin, zmax, dz);
					const Number_t abs_d = pat_wave_3d_absortion(bz+dz/2, zmin, zmax, pml_strength, pml_width);
					const Number_t dir_d = pat_wave_3d_directional(idt, abs_d);
					const Number_t upd_d = pat_wave_3d_pre_update(idt, abs_d, dir_d);
					const Number_t div_d = pat_wave_3d_pre_divergence(density, c, dir_d, dz);

					diver += div_d*(cache[tx][ty][3]-local_z);
					update = (upd_x + upd_y + upd_d)/3;
				}
				
				const Number_t local_new = update*cache[tx][ty][0] + diver;
				u[idx+0] = local_new;

				Number_t sint = sin(frequency*t);
				if(steps > 0){
					Number_t pamp = amplitude[base];
					
					Number_t kk = pamp*pamp*2;
					Number_t amp = sqrt((kk+((local_new*local_new - kk)/steps))/2);
					amplitude[base] = amp;
					if(abs(sint) > 0.2f && abs(sint) < 0.8f){
						const Number_t oldPha = phase[base]; 
						const Number_t pha = pat_phase(cache[tx][ty][0], local_new, t-dt, t, frequency);
						phase[base] = interpolate_angles(oldPha, pha, 0.9);
					}
				}
			}
		}
	}
}

__global__ void cuda_pat_wave_3d_listen_kernel(const Number_t * __restrict__ u,
											   Number_t * __restrict__ output,
										       const int num_listening,
											   const Number_t * __restrict__ listening_positions,
											   const int nx,
											   const int ny,
											   const int nz){

	const Number_t dx = kernel_constants[2];
	const Number_t dy = kernel_constants[3];
	const Number_t xmin = kernel_constants[4];
	const Number_t xmax = kernel_constants[5];
	const Number_t ymin = kernel_constants[6];
	const Number_t ymax = kernel_constants[7];
	const Number_t zmin = kernel_constants[8];
	const Number_t zmax = kernel_constants[9];
	const Number_t dz = kernel_constants[17];

	int ii = blockIdx.x*blockDim.x + threadIdx.x;

	if(ii < num_listening){
		Number_t data[2][2][2];
		Number_t x = listening_positions[3*ii];
		Number_t y = listening_positions[3*ii + 1];
		Number_t z = listening_positions[3*ii + 2];

		const int i = (int) floor((x - dx/2 - xmin)/dx);
		const int j = (int) floor((y - dy/2 - ymin)/dy);
		const int k = (int) floor((z - dz/2 - zmin)/dz);

		x = (x - w_get_pos(i, nx, xmin, xmax, dx))/dx;
		y = (y - w_get_pos(j, ny, ymin, ymax, dy))/dy;
		z = (z - w_get_pos(k, nz, zmin, zmax, dz))/dz;
	
		data[0][0][0] = u[4*(i + nx*(j + ny*k))];
		data[0][0][1] = u[4*((i+1) + nx*(j + ny*k))];
		data[0][1][0] = u[4*(i + nx*((j+1) + ny*k))];
		data[0][1][1] = u[4*((i+1) + nx*((j+1) + ny*k))];
		data[1][0][0] = u[4*(i + nx*(j + ny*(k+1)))];
		data[1][0][1] = u[4*((i+1) + nx*(j + ny*(k+1)))];
		data[1][1][0] = u[4*(i + nx*((j+1) + ny*(k+1)))];
		data[1][1][1] = u[4*((i+1) + nx*((j+1) + ny*(k+1)))];

		data[0][0][0] = (1-x)*data[0][0][0] + x*data[0][0][1];
		data[0][1][0] = (1-x)*data[0][1][0] + x*data[0][1][1];
		data[1][0][0] = (1-x)*data[1][0][0] + x*data[1][0][1];
		data[1][1][0] = (1-x)*data[1][1][0] + x*data[1][1][1];

		data[0][0][0] = (1-y)*data[0][0][0] + x*data[0][1][0];
		data[1][0][0] = (1-y)*data[1][0][0] + x*data[1][1][0];

		data[0][0][0] = (1-z)*data[0][0][0] + z*data[1][0][0];

		output[ii] = data[0][0][0];
	}
}

__device__ void cuda_pat_compute_multipole_term(Number_t * __restrict__ multipole,
												int base,
											    const Number_t * __restrict__ amplitude,
											   	const Number_t * __restrict__ phase,
											   	const int num_multipole,
											   	int nx, int ny, int nz,
											   	Number_t x, Number_t y, Number_t z,
											   	Number_t * __restrict__ legendre_base,
											   	Number_t p[2],
											   	Number_t dp_dn[2]){
	
	const Number_t xmin = kernel_constants[4];
	const Number_t xmax = kernel_constants[5];
	const Number_t ymin = kernel_constants[6];
	const Number_t ymax = kernel_constants[7];
	const Number_t zmin = kernel_constants[8];
	const Number_t zmax = kernel_constants[9];
	const Number_t c = kernel_constants[0];
	const Number_t frequency = kernel_constants[16];
	const Number_t ko = frequency/c;
	const Number_t dx = kernel_constants[2];
	const Number_t dy = kernel_constants[3];
	const Number_t dz = kernel_constants[17];
	const Number_t cx = x-kernel_constants[13];
	const Number_t cy = y-kernel_constants[14];
	const Number_t cz = z-kernel_constants[15];
	const Number_t r = sqrt(cx*cx + cy*cy + cz*cz);
	const Number_t theta = acos(cz / r);
    const Number_t phi  = atan2(cy, cx);
    const Number_t PI = acos(-1.0f);
    const Number_t area = (abs(dx*dy*cz) + abs(dx*cy*dz) + abs(cx*dy*dz))/r;

    const int idt = threadIdx.x;

    const int i = (int) floor((x - dx/2 - xmin)/dx);
	const int j = (int) floor((y - dy/2 - ymin)/dy);
	const int kk = (int) floor((z - dz/2 - zmin)/dz);

	//Compute p and dp_dn in the shared memory

	__syncthreads();

	if(idt == 0){
		Number_t t[3][2];
		Number_t ph = phase[i + nx*(j + ny*kk)];
		p[0] = cos(ph); p[1] = sin(ph);
		ph = phase[i+1 + nx*(j + ny*kk)];
		t[0][0] = cos(ph); t[0][1] = sin(ph);
		ph = phase[i + nx*(j+1 + ny*kk)];
		t[1][0] = cos(ph); t[1][1] = sin(ph);
		ph = phase[i + nx*(j + ny*(kk+1))];
		t[2][0] = cos(ph); t[2][1] = sin(ph);

		Number_t amp = amplitude[i + nx*(j + ny*kk)];
		p[0] *= amp; p[1] *= amp;
		amp = amplitude[i+1 + nx*(j + ny*kk)];
		t[0][0] *= amp; t[0][1] *= amp;
		amp = amplitude[i + nx*(j+1 + ny*kk)];
		t[1][0] *= amp; t[1][1] *= amp;
		amp = amplitude[i + nx*(j + ny*(kk+1))];
		t[2][0] *= amp; t[2][1] *= amp;

		dp_dn[0] = ((t[0][0]-p[0])*cx/dx + (t[1][0]-p[0])*cy/dy  + (t[2][0]-p[0])*cz/dz)/r;
		dp_dn[1] = ((t[0][1]-p[1])*cx/dx + (t[1][1]-p[1])*cy/dy  + (t[2][1]-p[1])*cz/dz)/r;
	}

	__syncthreads();
	//Compute the P_mm in the shared memory
	const Number_t xx = cos(theta);
	const Number_t somx2 = sqrt((1-xx)*(1+xx));

	if(idt == 0){
		legendre_base[0] = 0;
		legendre_base[1] = sqrt(1/(4*PI));

		for(int m = 0; m <= num_multipole; m++){
			legendre_base[m+2] = -legendre_base[m+1]*sqrt((2*m+3.0f)/(2*m+2.0f))*somx2;
		}
	}
	__syncthreads();

	const int m = idt;
	Number_t legendre[3][3];

	//[m][l]
	legendre[0][0] = 0;
	legendre[0][1] = legendre_base[m];
	legendre[0][2] = legendre_base[m]*sqrt(2*m+1.0f)*xx;

	legendre[1][0] = legendre[1][1] = 0;
	legendre[1][2] = legendre_base[m+1];

	legendre[2][0] = legendre[2][1] = legendre[2][2] = 0;
	
	for(int l = m; l <= num_multipole; l++){
		legendre[0][0] = legendre[0][1]; legendre[0][1] = legendre[0][2];
		legendre[1][0] = legendre[1][1]; legendre[1][1] = legendre[1][2];
		legendre[2][0] = legendre[2][1]; legendre[2][1] = legendre[2][2];

		if(l == m){
			legendre[2][2] = legendre_base[m+2];
			legendre[1][2] = legendre[1][1]*sqrt(2*m+3.0f)*xx;
			legendre[0][2] = (m > 0 ? (xx*sqrt(((2*l+1.0f)*(2*l+3.0f))/((l-m+2.0f)*(l+m+0.0f)))*legendre[0][1] - sqrt(((2*l+3.0f)*(l-m+1.0f)*(l+m-1.0f))/((2*l-1.0f)*(l-m+2.0f)*(l+m+0.0f)))*legendre[0][0]) : 0.0f);
		} else if(l == m+1){
			legendre[2][2] = legendre[2][1]*sqrt(2*m+5.0f)*xx;
			legendre[1][2] = xx*sqrt(((2*l+1.0f)*(2*l+3.0f))/((l-m+1.0f)*(l+m+1.0f)))*legendre[1][1] - sqrt(((2*l+3.0f)*(l-m)*(l+m))/((2*l-1.0f)*(l-m+1.0f)*(l+m+1.0f)))*legendre[1][0];
			legendre[0][2] = xx*sqrt(((2*l+1.0f)*(2*l+3.0f))/((l-m+2.0f)*(l+m+0.0f)))*legendre[0][1] - sqrt(((2*l+3.0f)*(l-m+1.0f)*(l+m-1.0f))/((2*l-1.0f)*(l-m+2.0f)*(l+m+0.0f)))*legendre[0][0];
		} else{
			legendre[2][2] = xx*sqrt(((2*l+1.0f)*(2*l+3.0f))/((l-m+0.0f)*(l+m+2.0f)))*legendre[2][1] - sqrt(((2*l+3.0f)*(l-m-1.0f)*(l+m+1.0f))/((2*l-1.0f)*(l-m+0.0f)*(l+m+2.0f)))*legendre[2][0];
			legendre[1][2] = xx*sqrt(((2*l+1.0f)*(2*l+3.0f))/((l-m+1.0f)*(l+m+1.0f)))*legendre[1][1] - sqrt(((2*l+3.0f)*(l-m+0.0f)*(l+m+0.0f))/((2*l-1.0f)*(l-m+1.0f)*(l+m+1.0f)))*legendre[1][0];
			legendre[0][2] = xx*sqrt(((2*l+1.0f)*(2*l+3.0f))/((l-m+2.0f)*(l+m+0.0f)))*legendre[0][1] - sqrt(((2*l+3.0f)*(l-m+1.0f)*(l+m-1.0f))/((2*l-1.0f)*(l-m+2.0f)*(l+m+0.0f)))*legendre[0][0];
		}

		//Compute multipole l, +m -> Use R l, -m
		{
			Number_t contrib[2] = {0, 0};
			const int mm = -m;
			Number_t R[2] = {0, 0};
			regular_basis(mm, l, phi, 1, 1, R, legendre);
			complexMultAdd(R, dp_dn, contrib);

			Number_t dR_dn[2] = {0, 0};
			
			//Calculating dR_dN
			{
				Number_t A[2] = {0, 0};
				Number_t temp[2] = {0, 0};
				regular_basis(mm+1, l+1, phi, 2, 2, temp, legendre); complexScaleAdd(Bmn(-mm-1, l+1), temp, A);
				regular_basis(mm+1, l-1, phi, 2, 0, temp, legendre); complexScaleAdd(-Bmn(mm, l), temp, A);

				Number_t B[2] = {0, 0};
				regular_basis(mm-1, l+1, phi, 0, 2, temp, legendre); complexScaleAdd(Bmn(mm-1, l+1), temp, B);
				regular_basis(mm-1, l-1, phi, 0, 0, temp, legendre); complexScaleAdd(-Bmn(-mm, l), temp, B);

				//rx
				temp[0] = temp[1] = 0;
				complexScaleAdd(0.5f, A, temp); complexScaleAdd(0.5f, B, temp);
				complexScaleAdd(cx/r, temp, dR_dn);

				//ry
				temp[0] = (A[1]-B[1])*0.5f; temp[1] = (B[0] - A[0])*0.5f;
				complexScaleAdd(cy/r, temp, dR_dn);

				//Use A as rz
				A[0] = 0; A[1] = 0;
				regular_basis(mm, l-1, phi, 1, 0, temp, legendre); complexScaleAdd(Amn(mm, l-1), temp, A);
				regular_basis(mm, l+1, phi, 1, 2, temp, legendre); complexScaleAdd(-Amn(mm, l), temp, A);
				complexScaleAdd(cz/r, A, dR_dn);
				
				dR_dn[0] *= -ko; dR_dn[1] *= -ko;
			}
			complexMultAdd(p, dR_dn, contrib);
			contrib[0] *= area; contrib[1] *= area;
			// regular_basis(mm, l, phi, 1, 1, contrib, legendre);
			//contrib[0] = dR_dn[0]; contrib[1] = dR_dn[1]; 
			//Sum to the correspondent multipole coefficient (l*(l+1)+m)
			multipole[base + 2*(l*(l+1) + m)] += contrib[0];
			multipole[base + 2*(l*(l+1) + m)+1] += contrib[1];
		}

		//Compute multipole l, -m -> Use R l, m
		if(m != 0){
			Number_t contrib[2] = {0, 0};
			const int mm = m;
			Number_t R[2] = {0, 0};
			regular_basis(mm, l, phi, 1, 1, R, legendre);
			complexMultAdd(R, dp_dn, contrib);

			Number_t dR_dn[2] = {0, 0};
			
			//Calculating dR_dN
			{
				Number_t A[2] = {0, 0};
				Number_t temp[2] = {0, 0};
				regular_basis(mm+1, l+1, phi, 2, 2, temp, legendre); complexScaleAdd(Bmn(-mm-1, l+1), temp, A);
				regular_basis(mm+1, l-1, phi, 2, 0, temp, legendre); complexScaleAdd(-Bmn(mm, l), temp, A);

				Number_t B[2] = {0, 0};
				regular_basis(mm-1, l+1, phi, 0, 2, temp, legendre); complexScaleAdd(Bmn(mm-1, l+1), temp, B);
				regular_basis(mm-1, l-1, phi, 0, 0, temp, legendre); complexScaleAdd(-Bmn(-mm, l), temp, B);

				//rx
				temp[0] = temp[1] = 0;
				complexScaleAdd(0.5f, A, temp); complexScaleAdd(0.5f, B, temp);
				complexScaleAdd(cx/r, temp, dR_dn);

				temp[0] = (A[1]-B[1])*0.5f; temp[1] = (B[0] - A[0])*0.5f;
				complexScaleAdd(cy/r, temp, dR_dn);

				//Use A as rz
				A[0] = 0; A[1] = 0;
				regular_basis(mm, l-1, phi, 1, 0, temp, legendre); complexScaleAdd(Amn(mm, l-1), temp, A);
				regular_basis(mm, l+1, phi, 1, 2, temp, legendre); complexScaleAdd(-Amn(mm, l), temp, A);

				complexScaleAdd(cz/r, A, dR_dn);
				dR_dn[0] *= -ko; dR_dn[1] *= -ko;
			}
			complexMultAdd(p, dR_dn, contrib);
			contrib[0] *= area; contrib[1] *= area;
			// regular_basis(mm, l, phi, 1, 1, contrib, legendre);
			// contrib[0] = dR_dn[0]; contrib[1] = dR_dn[1]; 
			//Sum to the correspondent multipole coefficient (l*(l+1)-m)
			multipole[base + 2*(l*(l+1) - m)] += contrib[0];
			multipole[base + 2*(l*(l+1) - m)+1] += contrib[1];
		}
	}
}

__global__ void cuda_pat_wave_3d_compute_multipoles(Number_t * __restrict__ integral_multipole,
													const Number_t * __restrict__ amplitude,
											   		const Number_t * __restrict__ phase,
											        const int num_multipole,
												    const int nx,
												    const int ny,
												    const int nz,
												    const Number_t radius){

	//Decide which points I am going to compute
	const Number_t xmin = kernel_constants[4];
	const Number_t xmax = kernel_constants[5];
	const Number_t ymin = kernel_constants[6];
	const Number_t ymax = kernel_constants[7];
	const Number_t zmin = kernel_constants[8];
	const Number_t zmax = kernel_constants[9];
	const Number_t cx = kernel_constants[13];
	const Number_t cy = kernel_constants[14];
	const Number_t cz = kernel_constants[15];
	const Number_t dx = kernel_constants[2];
	const Number_t dy = kernel_constants[3];
	const Number_t dz = kernel_constants[17];

	__shared__ Number_t legendre_base[200];
	__shared__ Number_t dp_dn[2];
	__shared__ Number_t p[2];

	const int i = blockIdx.x;
	const int idt = threadIdx.x;
	if(i < nx){
		const int nn = (num_multipole+1)*(num_multipole+1);
		const int base = i*2*nn;
		if(idt == 0){
			for(int it = 0; it < nn; it++){
				integral_multipole[base+2*it] = 0;
				integral_multipole[base+2*it+1] = 0;
			}
		}
		__syncthreads();

		const Number_t bx = w_get_pos(i, nx, xmin, xmax, dx);
		for(int k = 0; k < nz; k++){
			const Number_t bz = w_get_pos(k, nz, zmin, zmax, dz);
			const Number_t sqdist = radius*radius - (bx-cx)*(bx-cx) - (bz-cz)*(bz-cz);
			if(sqdist > 0){
				const Number_t ya = cy + sqrt(sqdist);
				const Number_t yb = cy - sqrt(sqdist);

				cuda_pat_compute_multipole_term(integral_multipole,
												base,
											    amplitude,
											   	phase,
											   	num_multipole,
											   	nx, ny, nz,
											   	bx, ya, bz,
											   	legendre_base,
											   	p, dp_dn);

				cuda_pat_compute_multipole_term(integral_multipole,
												base,
											    amplitude,
											   	phase,
											   	num_multipole,
											   	nx, ny, nz,
											   	bx, yb, bz,
											   	legendre_base,
											   	p, dp_dn);
			} else if(sqdist == 0){
				cuda_pat_compute_multipole_term(integral_multipole,
												base,
											    amplitude,
											   	phase,
											   	num_multipole,
											   	nx, ny, nz,
											   	bx, cy, bz,
											   	legendre_base,
											   	p, dp_dn);
			}
		}
	}
}

__global__ void cuda_pat_wave_3d_reduce_multipoles(Number_t * __restrict__ multipole,
											   	   const Number_t * __restrict__ integral_multipole,
											       const int num_multipole,
												   const int nx){
	const int ii = blockIdx.x*blockDim.x + threadIdx.x;
	const int nn = (num_multipole+1)*(num_multipole+1);
	const Number_t c = kernel_constants[0];
	const Number_t frequency = kernel_constants[16];
	const Number_t k = frequency/c;

	if(ii < nn){
		Number_t coef[2] = {0, 0};
		int idx = 0;
		for(int i = 0; i < nx; i++){
			coef[0] += integral_multipole[idx + 2*ii];
			coef[1] += integral_multipole[idx + 2*ii + 1];
			idx += 2*nn;
		}
		// multipole[2*ii] = coef[0];
		// multipole[2*ii + 1] = coef[1];
		multipole[2*ii] = -k*coef[1];
		multipole[2*ii + 1] = k*coef[0];
	}
}

#undef RRADX
#undef RRADY	