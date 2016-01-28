#include <cstdio>
#include <cmath>
#include <cstdlib>

#ifdef USE_CUDA
	#include "cuda/cuda_PAN_wave_3d.h"
#else
	#include "cpu/cpu_pml_wave_2d.h"
#endif


Number_t zeros(const Number_t x, const Number_t y, const Number_t z){
	return 0;
}

bool sphere(const Number_t x, const Number_t y, const Number_t z){
	Number_t d = sqrt(  (x-0.5)*(x-0.5)
					  + (y-0.5)*(y-0.5)
					  + (z-0.5)*(z-0.5));
	return (d < 0.1);
}

Number_t sphereGradient(const Number_t x, const Number_t y, const Number_t z, int dim){
	Number_t d = sqrt(  (x-0.5)*(x-0.5)
					  + (y-0.5)*(y-0.5)
					  + (z-0.5)*(z-0.5));

	if(dim == 0){
		return (x-0.5)/d;
	} else if(dim == 1){
		return (y-0.5)/d;
	} else{
		return (z-0.5)/d;
	}
}

void writeToFile(FILE * fp, Number_t * u, int nx, int ny){
	for(int i = 0; i < ny; i++){
		for(int j =0; j < nx; j++){
			fprintf(fp, "%.5f ", u[3*(j + i*nx)]);
		}
		fprintf(fp, "\n");
	}
}

int main(int argc, char** argv){
	#ifdef USE_CUDA
		Cuda_PAN_Wave_3d_t wave;
	#else
		Cpu_PML_Wave_2d_t wave;
	#endif

	int nx = 176;
	Number_t cellsize = 1.0/nx;
	char filename[1024];
	printf("%.3lf/n", 6*sizeof(Number_t)*(nx)*(nx)*(nx)/(1024.0*1024.0*1024.0));
	int nsteps = 500;
	Number_t c = 0.34029;
	Number_t dt = 0.6/(nx*c);
	Number_t * u = NULL;
	int lis = 20000;
	Number_t * listening = (Number_t *) malloc(3*lis*sizeof(Number_t));
	for(int i = 0; i < lis; i++){
		Number_t x = (Number_t)rand()/RAND_MAX; x = 0.7*x + 0.15;
		Number_t y = (Number_t)rand()/RAND_MAX; y = 0.7*y + 0.15;
		Number_t z = (Number_t)rand()/RAND_MAX; z = 0.7*z + 0.15;

		listening[3*i] = x;
		listening[3*i + 1] = y;
		listening[3*i + 2] = z;
	}
	wave = wave_sim_init(0, 0, 0,
						 1, 1, 1,
						 c, 0.6/(nx*c),
						 cellsize,
						 lis,
						 listening,
						 &zeros,
						 &sphere,
						 0.5, 0.5, 0.5,
						 &sphereGradient,
						 0.1,
						 100,
						 dt*10);


	for(int step = 0; step < nsteps; step++){
		//u = wave_sim_get_u(wave);
		printf("Frame %d %lf\n", step, 0.6/(nx*c));
		sprintf(filename, "frames/frame%d", step);
		// FILE *fp = fopen(filename, "w+");
		// writeToFile(fp, u, nx, ny);
		// fclose(fp);
		wave_sim_step(wave);
		wave_listen(wave, 0);
		wave_listen(wave, 1);
		wave_listen(wave, 2);
		wave_listen(wave, 3);
		wave_listen(wave, 4);
		wave_listen(wave, 5);
	}

	// 	u = wave_sim_get_u(wave, 0);
	// for(int i = 0; i < ny+2; i++){
	// 	for(int j = 0; j < nx+2; j++){
	// 		printf("%.4f\t\t", u[j + i*(nx+2)]);
	// 	}
	// 	printf("\n");
	// }
	// wave_sim_step(wave);
	// u = wave_sim_get_u(wave, 0);
	// for(int i = 0; i < ny+2; i++){
	// 	for(int j = 0; j < nx+2; j++){
	// 		printf("%.3f\t\t", u[j + i*(nx+2)]);
	// 	}
	// 	printf("\n");
	// }
}