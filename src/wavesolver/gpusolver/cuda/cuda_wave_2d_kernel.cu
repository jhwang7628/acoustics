__global__ void cuda_wave_2d_kernel(Number_t * __restrict__ unew,
									const Number_t * __restrict__ u,
									const Number_t * __restrict__ uold,
									const int nx, const int ny,
									const Number_t c, const Number_t dt,
									const Number_t dx, const Number_t dy){
	#define BDIMX 16
	#define BDIMY 16

	__shared__ Number_t cache[BDIMX + 2][BDIMY + 2];

	int j = blockIdx.x*blockDim.x + threadIdx.x + 1;
	int i = blockIdx.y*blockDim.y + threadIdx.y + 1;

	if(j <= nx && i <= ny){

		int idx = j + (nx+2)*i;

		 int bdy = min(BDIMY, ny-i+1);
		 int bdx = min(BDIMX, nx-j+1);
		int tx = threadIdx.x + 1;
		int ty = threadIdx.y + 1;

		Number_t oldU = uold[idx];
		Number_t unow = u[idx];
		__syncthreads();
		//cache[x][y]
		if(threadIdx.y < 1){
			cache[tx][threadIdx.y] = u[j + (nx+2)*(i-1)];
			cache[tx][threadIdx.y+bdy+1] = u[j + (nx+2)*(i+bdy)];
		}
		if(threadIdx.x < 1){
			cache[threadIdx.x][ty] = u[(j-1) + (nx+2)*i];
			cache[threadIdx.x+bdx+1][ty] = u[(j+bdx) + (nx+2)*i];
		}
		cache[tx][ty] = unow;
		__syncthreads();

		Number_t tauX = c*dt/dx;
		Number_t tauY = c*dt/dy;

		Number_t tauX2 = tauX*tauX;
		Number_t tauY2 = tauY*tauY;

		Number_t newU = 2*cache[tx][ty] - oldU
					  + tauY2*(cache[tx][ty+1] - 2*cache[tx][ty] + cache[tx][ty-1])
					  + tauX2*(cache[tx+1][ty] - 2*cache[tx][ty] + cache[tx-1][ty])
					  ;

		unew[idx] = newU;
	}
}