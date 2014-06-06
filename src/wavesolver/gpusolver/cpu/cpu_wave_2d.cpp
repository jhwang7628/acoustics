#include "cpu_wave_2d.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cstdio>

struct Cpu_Wave_2d_sim_data_t {
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
};



Cpu_Wave_2d_t wave_sim_init(Number_t xmin, Number_t ymin,
						Number_t xmax, Number_t ymax,
						Number_t c, Number_t dt, 
						int nx, int ny,
						Number_t (*init_function)(Number_t, Number_t, void *),
						void * ctx){
	Cpu_Wave_2d_t wave = (Cpu_Wave_2d_t) malloc(sizeof(Cpu_Wave_2d_sim_data_t));
	
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

	wave->ubuf = (Number_t *) malloc(3*(wave->nx+2)*(wave->ny+2)*sizeof(Number_t));

	Number_t * u = wave_sim_get_u(wave, 0);
	for(int i = 1; i <= ny; i++){
		Number_t y = wave_sim_get_y(wave, i);
		for(int j = 1; j <= nx; j++){
			Number_t x = wave_sim_get_x(wave, j);
			u[j + i*(nx+2)] = init_function(x, y, ctx);
		}
	}

	wave_sim_apply_boundary(wave);

	Number_t * u_old = wave_sim_get_u(wave, -1);
	memcpy(u_old, u, (wave->nx+2)*(wave->ny+2)*sizeof(Number_t));

	return wave;
}

void wave_sim_free(Cpu_Wave_2d_t wave){
	free(wave->ubuf);
	free(wave);
}

Number_t wave_sim_get_x(Cpu_Wave_2d_t wave, int j){
	return ((j-1)*wave->xmax + (wave->nx - j)*wave->xmin)/(wave->nx-1);
}

Number_t wave_sim_get_y(Cpu_Wave_2d_t wave, int i){
	return ((i-1)*wave->ymax + (wave->ny - i)*wave->ymin)/(wave->ny-1);
}

Number_t * wave_sim_get_u(Cpu_Wave_2d_t wave, int offset){
	return wave->ubuf + ((wave->step + offset + 3)%3)*(wave->nx+2)*(wave->ny+2);
}

void wave_sim_step(Cpu_Wave_2d_t wave){
	clock_t t;

	t = clock();
	
		Number_t * uold = wave_sim_get_u(wave, -1);
		Number_t * u 	= wave_sim_get_u(wave,  0);
		Number_t * unew = wave_sim_get_u(wave,  1);

		int nx = wave->nx;
		int ny = wave->ny;

		Number_t tauX = wave->c*wave->dt/wave->dx;
		Number_t tauY = wave->c*wave->dt/wave->dy;

		Number_t tauX2 = tauX*tauX;
		Number_t tauY2 = tauY*tauY;

		Number_t cache[3][3];

		for(int i = 1; i <= ny; i++){
			cache[0][0] = u[(1-1)  + (i-1)*(nx+2)];
			cache[0][1] = u[(1)    + (i-1)*(nx+2)];
			cache[1][0] = u[(1-1)  + (i)*(nx+2)];
			cache[1][1] = u[(1)    + (i)*(nx+2)];
			cache[2][0] = u[(1-1)  + (i+1)*(nx+2)];
			cache[2][1] = u[(1)    + (i+1)*(nx+2)];
			
			for(int j = 1; j <= nx; j++){
				cache[0][2] = u[(j+1)  + (i-1)*(nx+2)];
				cache[1][2] = u[(j+1)  + (i)*(nx+2)];
				cache[2][2] = u[(j+1) + (i+1)*(nx+2)];

				//compute
				Number_t oldU = uold[j + i*(nx+2)];

				Number_t newU = 2*cache[1][1] - oldU
							  + (tauX2)*(cache[1][2] - 2*cache[1][1] + cache[1][0])
							  + (tauY2)*(cache[2][1] - 2*cache[1][1] + cache[0][1]);

				unew[j + i*(nx+2)] = newU;

				cache[0][0] = cache[0][1]; cache[0][1] = cache[0][2];
				cache[1][0] = cache[1][1]; cache[1][1] = cache[1][2];
				cache[2][0] = cache[2][1]; cache[2][1] = cache[2][2];
			}
		}
	t = clock()-t;
	printf("Step Time: %fms\n", 1000*((float)t)/CLOCKS_PER_SEC);
	
		wave->step++;
		wave->t += wave->dt;
	t = clock();
		wave_sim_apply_boundary(wave);
	t = clock()-t;
	printf("Boundary Time: %fms\n", 1000*((float)t)/CLOCKS_PER_SEC);
	
}

void wave_sim_apply_boundary(Cpu_Wave_2d_t wave){
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