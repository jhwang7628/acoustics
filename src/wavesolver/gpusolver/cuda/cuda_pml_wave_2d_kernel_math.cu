__device__ __forceinline__ Number_t pml_wave_2d_absortion(const Number_t pos, const Number_t dmin, const Number_t dmax, const Number_t strength, const Number_t width){
	Number_t d;
	if(pos-dmin < width){
		d = (width-(pos-dmin))/width;
		return strength*d*d;
	} else if(dmax-pos < width){
		d = (width-(dmax-pos))/width;
		return strength*d*d;
	} else{
		return 0;
	}
}

__device__ __forceinline__ Number_t pml_wave_2d_vel_update(const Number_t idt, const Number_t absortion){
	return (idt - absortion/2.0)/(idt + absortion/2.0);
}

__device__ __forceinline__ Number_t pml_wave_2d_directional(const Number_t idt, const Number_t absortion){
	return (idt + absortion/2.0);
}

__device__ __forceinline__ Number_t pml_wave_2d_pre_update(const Number_t idt, const Number_t absortion, const Number_t directional){
	return (idt - absortion/2.0)/directional;
}

__device__ __forceinline__ Number_t pml_wave_2d_pre_divergence(const Number_t density, const Number_t c, const Number_t directional, const Number_t cellSize){
	return (-1.0*density*(c/cellSize)*(c/directional));
}

__device__ __forceinline__ Number_t pml_wave_2d_gradient(const Number_t idt, const Number_t absortion, const Number_t cellSize, const Number_t density){
	return -1.0/(density*cellSize*(idt + absortion/2.0));
}