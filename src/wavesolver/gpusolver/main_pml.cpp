#include <cstdio>
#include <cmath>

#ifdef USE_CUDA
	#include "cuda/cuda_pml_wave_2d.h"
#else
	#include "cpu/cpu_pml_wave_2d.h"
#endif


Number_t gaussian(Number_t x, Number_t y, void * ctx){
	Number_t stddev = 0.01;
	Number_t mean = 0.2;
	Number_t var2 = stddev*stddev*2;
	Number_t term = sqrt((x-mean)*(x-mean) + (y-mean)*(y-mean));
	// Number_t term = x-mean;
	return stddev*exp(-term*term/var2)/sqrt(acos(-1)*var2);
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
		Cuda_PML_Wave_2d_t wave;
	#else
		Cpu_PML_Wave_2d_t wave;
	#endif

	int nx = 256;
	char filename[1024];
	int ny = nx;
	printf("%.3lf/n", 3*sizeof(Number_t)*(nx)*(nx)/(1024.0*1024.0*1024.0));
	int nsteps = 200;
	Number_t c = 0.34029;
	Number_t * u = NULL;
	wave = wave_sim_init(0, 0, 1, 1,
						c, 0.6/(nx*c),
						nx, ny,
						gaussian,
						NULL,
						// 0.1,
						// 0.0
						0.1,
						100.0
						);


	for(int step = 0; step < nsteps; step++){
		u = wave_sim_get_u(wave);
		printf("Frame %d\n", step);
		sprintf(filename, "frames/frame%d", step);
		FILE *fp = fopen(filename, "w+");
		writeToFile(fp, u, nx, ny);
		fclose(fp);
		wave_sim_step(wave);
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