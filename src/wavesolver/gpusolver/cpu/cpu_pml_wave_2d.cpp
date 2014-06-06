#include "cpu_pml_wave_2d.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cstdio>

struct Cpu_PML_Wave_2d_sim_data_t {
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
};

#include "cpu_pml_wave_2d_math.cpp"

Number_t * wave_sim_get_u(Cpu_PML_Wave_2d_t wave){
	return wave->ubuf;
}

Cpu_PML_Wave_2d_t wave_sim_init(Number_t xmin, Number_t ymin,
						Number_t xmax, Number_t ymax,
						Number_t c, Number_t dt, 
						int nx, int ny,
						Number_t (*init_function)(Number_t, Number_t, void *),
						void * ctx,
						Number_t pml_width,
						Number_t pml_strength){
	Cpu_PML_Wave_2d_t wave = (Cpu_PML_Wave_2d_t) malloc(sizeof(Cpu_PML_Wave_2d_sim_data_t));
	
	wave->xmin = xmin;
	wave->ymin = ymin;
	wave->xmax = xmax;
	wave->ymax = ymax;

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

	wave->ubuf = (Number_t *) malloc(3*(wave->nx)*(wave->ny)*sizeof(Number_t));

	//Set the pressures
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
		if(i == ny-1){
			for(int j = 0; j < nx; j++){
				if(j == nx-1){
					u[3*(j + i*nx) + 1] = 0;
					u[3*(j + i*nx) + 2] = 0;
				} else{
					u[3*(j + i*nx) + 1] = -(u[3*(j+1 + i*nx)] - u[3*(j + i*nx)])*wave->dt/wave->dx;
					u[3*(j + i*nx) + 2] = 0;
				}
			}
		} else{
			for(int j = 0; j < nx; j++){
				if(j == nx-1){
					u[3*(j + i*nx) + 1] = 0;
					u[3*(j + i*nx) + 2] = -(u[3*(j + (i+1)*nx)] - u[3*(j + i*nx)])*wave->dt/wave->dy;
				} else{
					u[3*(j + i*nx) + 1] = -(u[3*(j+1 + i*nx)] - u[3*(j + i*nx)])*wave->dt/wave->dx;
					u[3*(j + i*nx) + 2] = -(u[3*(j + (i+1)*nx)] - u[3*(j + i*nx)])*wave->dt/wave->dy;
				}
			}
		}
	}


	wave_sim_apply_boundary(wave);

	return wave;
}

void wave_sim_free(Cpu_PML_Wave_2d_t wave){
	free(wave->ubuf);
	free(wave);
}

Number_t wave_sim_get_x(Cpu_PML_Wave_2d_t wave, int j){
	return ((j*wave->xmax + (wave->nx - j)*wave->xmin)/wave->nx) + wave->dx/2;
}

Number_t wave_sim_get_y(Cpu_PML_Wave_2d_t wave, int i){
	return ((i*wave->ymax + (wave->ny - i)*wave->ymin)/wave->ny) + wave->dy/2;
}

void wave_sim_step(Cpu_PML_Wave_2d_t wave){
	clock_t t;

	t = clock();
		//Update V
		Number_t * u = wave->ubuf;


		int nx = wave->nx;
		int ny = wave->ny;

		Number_t cache[2][2];
		for(int i = 0; i < ny; i++){
			cache[0][0] = u[3*(0 + i*nx)];
			if(i != ny-1) cache[0][1] = u[3*(0 + (i+1)*nx)]; 
			Number_t by = wave_sim_get_y(wave, i);
			for(int j = 0; j < nx; j++){
				Number_t oldVx = u[3*(j + i*nx)+1];
				Number_t oldVy = u[3*(j + i*nx)+2];
				cache[1][0] = u[3*((j+1) + i*nx)];
				if(j != nx-1) cache[1][1] = u[3*((j+1) + (i+1)*nx)]; 
				
				Number_t bx = wave_sim_get_x(wave, j);

				Number_t absortion;
				Number_t update;
				Number_t gradient;
				Number_t newVx = 0;
				Number_t newVy = 0;
				//X
				if(j != nx-1){
					absortion = pml_wave_2d_absortion(bx+wave->dx/2, wave->xmin, wave->xmax, wave->pml_strength, wave->pml_width);
					update = pml_wave_2d_vel_update(wave->dt, absortion);
					gradient = pml_wave_2d_gradient(wave->dt, absortion, wave->dx, wave->density);

					newVx = oldVx*update + gradient*(cache[1][0]-cache[0][0]);
				}

				//Y
				if(i != ny-1){
					absortion = pml_wave_2d_absortion(by+wave->dy/2, wave->ymin, wave->ymax, wave->pml_strength, wave->pml_width);
					update = pml_wave_2d_vel_update(wave->dt, absortion);
					gradient = pml_wave_2d_gradient(wave->dt, absortion, wave->dy, wave->density);

					newVy = oldVy*update + gradient*(cache[0][1]-cache[0][0]);
				}

				u[3*(j + i*nx)+1] = newVx;
				u[3*(j + i*nx)+2] = newVy;

				cache[0][0] = cache[1][0];
				cache[0][1] = cache[1][1];
			}
		}
	t = clock()-t;

	t = clock();
		//Update P
		//00 = down
		//10 = up
		//01 = left
		//11 = right
		Number_t abs_x;
		Number_t abs_y;
		Number_t dir_x;
		Number_t dir_y;
		Number_t upd_x;
		Number_t upd_y;
		Number_t div_x;
		Number_t div_y;
		for(int i = 0; i < ny; i++){
			Number_t by = wave_sim_get_y(wave, i);
			cache[0][1] = 0;
			cache[0][0] = 0;
			for(int j = 0; j < nx; j++){
				Number_t bx = wave_sim_get_x(wave, j);
				Number_t oldU = u[3*(j + i*nx)];
				cache[1][1] = u[3*(j + i*nx)+1];
				cache[1][0] = u[3*(j + i*nx)+2];
				if(i != 0) cache[0][0] = u[3*(j + (i-1)*nx)+2];

				abs_x = pml_wave_2d_absortion(bx+wave->dx/2, wave->xmin, wave->xmax, wave->pml_strength, wave->pml_width);
				abs_y = pml_wave_2d_absortion(by+wave->dy/2, wave->ymin, wave->ymax, wave->pml_strength, wave->pml_width);
				dir_x = pml_wave_2d_directional(wave->dt, abs_x);
				dir_y = pml_wave_2d_directional(wave->dt, abs_y);
				upd_x = pml_wave_2d_pre_update(wave->dt, abs_x, dir_x);
				upd_y = pml_wave_2d_pre_update(wave->dt, abs_y, dir_y);
				div_x = pml_wave_2d_pre_divergence(wave->density, wave->c, dir_x, wave->dx);
				div_y = pml_wave_2d_pre_divergence(wave->density, wave->c, dir_y, wave->dy);

				Number_t newU = oldU*(upd_x + upd_y)/2
							  + div_x*(cache[1][1]-cache[0][1])
							  + div_y*(cache[1][0]-cache[0][0]);
							  ;

				// newU = 2*oldU
				// + div_x*(cache[1][1]-cache[0][1])
				// + div_y*(cache[1][0]-cache[0][0]);
				// ;

				u[3*(j + i*nx)] = newU;

				cache[0][1] = cache[1][1];
			}
		}
	t = clock()-t;
	printf("Step Time: %fms\n", 1000*((float)t)/CLOCKS_PER_SEC);
	
		wave->t += wave->dt;
	t = clock();
		wave_sim_apply_boundary(wave);
	t = clock()-t;
	printf("Boundary Time: %fms\n", 1000*((float)t)/CLOCKS_PER_SEC);
	
}

void wave_sim_apply_boundary(Cpu_PML_Wave_2d_t wave){
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