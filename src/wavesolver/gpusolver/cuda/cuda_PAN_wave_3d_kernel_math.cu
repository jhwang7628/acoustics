__device__ __forceinline__ Number_t pan_wave_3d_absortion(const Number_t pos, const Number_t dmin, const Number_t dmax, const Number_t strength, const Number_t width){
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

__device__ __forceinline__ Number_t pan_wave_3d_vel_update(const Number_t idt, const Number_t absortion){
	return (idt - absortion/2.0)/(idt + absortion/2.0);
}

__device__ __forceinline__ Number_t pan_wave_3d_directional(const Number_t idt, const Number_t absortion){
	return (idt + absortion/2.0);
}

__device__ __forceinline__ Number_t pan_wave_3d_pre_update(const Number_t idt, const Number_t absortion, const Number_t directional){
	return (idt - absortion/2.0)/directional;
}

__device__ __forceinline__ Number_t pan_wave_3d_pre_divergence(const Number_t density, const Number_t c, const Number_t directional, const Number_t cellSize){
	return (-1.0*density*(c/cellSize)*(c/directional));
}

__device__ __forceinline__ Number_t pan_wave_3d_gradient(const Number_t idt, const Number_t absortion, const Number_t cellSize, const Number_t density){
	return -1.0/(density*cellSize*(idt + absortion/2.0));
}

__device__ Number_t PAN_Mitchelli(const Number_t t, const Number_t h){
	Number_t ret = 0;
	if(t <= 4*h){
		Number_t x = (t - 2*h)/h;
		x = x > 0? x : -x;
		
		if(x < 1.0){
			x = 1 - x;
			ret = -15*x;
			ret = x*(18 + ret);
			ret = x*(9 + ret);
			ret = 2 + ret;
		} else if(x < 2.0){
			x = 2-x;
			ret = 5*x;
			ret = x*x*(ret - 3);
		} else{
			ret = 0;
		}
		ret /= 18;
	}
	return ret;
}

__device__ Number_t PAN_boundary(const Number_t x, const Number_t y, const Number_t z, const int field, const int component){
	const Number_t cx = x-kernel_constants[13];
	const Number_t cy = y-kernel_constants[14];
	const Number_t cz = z-kernel_constants[15];

	Number_t retx=0, rety=0, retz=0;

	if(field == 0){
		retx = 1;
	} else if(field == 1){
		rety = 1;
	} else if(field == 2){
		retz = 1;
	} else if(field == 3){
		rety = cz;
		retz = -cy; 
	} else if(field == 4){
		retx = -cz;
		retz = cx;
	} else if(field == 5){
		retx = cy;
		rety = -cx;
	}

	if(component == 0){
		return retx;
	} else if(component == 1){
		return rety;
	} else{
		return retz;
	}
}