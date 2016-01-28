Number_t pml_wave_2d_absortion(Number_t pos, Number_t dmin, Number_t dmax, Number_t strength, Number_t width){
	Number_t h = pos-dmin;
	if(h < width){
		Number_t d = (width-h)/width;
		return strength*d*d;
	}
	h = dmax-pos;
	if(h < width){
		Number_t d = (width-h)/width;
		return strength*d*d;
	}
	return 0.0;
}

Number_t pml_wave_2d_vel_update(Number_t dt, Number_t absortion){
	return (1/dt - absortion/2.0)/(1/dt + absortion/2.0);
}

Number_t pml_wave_2d_directional(Number_t dt, Number_t absortion){
	return (1/dt + absortion/2.0);
}

Number_t pml_wave_2d_pre_update(Number_t dt, Number_t absortion, Number_t directional){
	return (1/dt - absortion/2.0)/directional;
}

Number_t pml_wave_2d_pre_divergence(Number_t density, Number_t c, Number_t directional, Number_t cellSize){
	return (-1.0*density*(c/cellSize)*(c/directional));
}

Number_t pml_wave_2d_gradient(Number_t dt, Number_t absortion, Number_t cellSize, Number_t density){
	return -1.0/(density*cellSize*(1/dt + absortion/2.0));
}