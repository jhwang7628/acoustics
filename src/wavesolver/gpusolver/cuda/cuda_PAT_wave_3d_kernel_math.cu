__device__ __forceinline__ Number_t pat_wave_3d_absortion(const Number_t pos, const Number_t dmin, const Number_t dmax, const Number_t strength, const Number_t width){
	Number_t d;
	Number_t pd = min(pos-dmin, dmax-pos);
	if(pd < width){
		d = (width-pd)/width;
		return strength*d*d;
	} else{
		return 0;
	}
}

__device__ __forceinline__ Number_t pat_wave_3d_vel_update(const Number_t idt, const Number_t absortion){
	return (idt - absortion/2.0f)/(idt + absortion/2.0f);
}

__device__ __forceinline__ Number_t pat_wave_3d_directional(const Number_t idt, const Number_t absortion){
	return (idt + absortion/2.0f);
}

__device__ __forceinline__ Number_t pat_wave_3d_pre_update(const Number_t idt, const Number_t absortion, const Number_t directional){
	return (idt - absortion/2.0f)/directional;
}

__device__ __forceinline__ Number_t pat_wave_3d_pre_divergence(const Number_t density, const Number_t c, const Number_t directional, const Number_t cellSize){
	return (-1.0f*density*(c/cellSize)*(c/directional));
}

__device__ __forceinline__ Number_t pat_wave_3d_gradient(const Number_t idt, const Number_t absortion, const Number_t cellSize, const Number_t density){
	return -1.0f/(density*cellSize*(idt + absortion/2.0f));
}

__device__ __forceinline__ Number_t pat_amplitude(const Number_t p1, const Number_t p2, const Number_t t1, const Number_t t2, const Number_t omega, const Number_t phase){
	Number_t cosi = sin(omega*t1 + phase);
	// if(cosi == 0.0f){
	// 	cosi = sin(omega*t2 + phase);
	// 	return abs(p2/cosi);
	// }
	return abs(p1/cosi);
}

__device__ __forceinline__ Number_t pat_phase(const Number_t p1, const Number_t p2, const Number_t t1, const Number_t t2, const Number_t omega){
	return atan2(p1*cos(omega*t2) - p2*cos(omega*t1), p2*sin(omega*t1) - p1*sin(omega*t2));
}

__device__ __forceinline__ Number_t Amn(int m, int n){
	int mm = abs(m);
    return n < mm ? 0 : sqrt(((mm+n+1.0f)*(n-mm+1.0f)) / ((2*n+1.0f)*(2*n+3.0f)));
}

__device__ __forceinline__ Number_t Bmn(int m, int n){
    return n < abs(m) ? 0 :
        (m >= 0 ? (sqrt(((n-m-1.0f)*(n-m)) / ((2*n-1.0f)*(2*n+1.0f))))
        		:-(sqrt(((n-m-1.0f)*(n-m)) / ((2*n-1.0f)*(2*n+1.0f))))
        		);
}

__device__ __forceinline__ Number_t bess(int n){
	return (n >= 0 ? bessel[n] : 0);
}

__device__ __forceinline__ void regular_basis(const int m, const int n, const Number_t phi, const int offm, const int offn, Number_t out[2], const Number_t legendre[3][3]){
	const Number_t sphi = sin(m*phi);
	const Number_t cphi = cos(m*phi);
	int kofm (m < 0 ? 2-offm : offm);
	const Number_t harm = ((n >= abs(m)) ? ((m&1) ? -(bess(n)*(legendre[kofm][offn]))
										     	  	   : (bess(n)*(legendre[kofm][offn]))
										     	       )
									     : 0
									     );

	out[0] = harm*cphi;
	out[1] = harm*sphi;
}

__device__ __forceinline__ void complexScaleAdd(const Number_t a, const Number_t b[2], Number_t out[2]){
	out[0] += a*b[0];
	out[1] += a*b[1];
}

__device__ __forceinline__ void complexMultAdd(const Number_t a[2], const Number_t b[2], Number_t out[2]){
	out[0] += a[0]*b[0]-a[1]*b[1];
	out[1] += a[0]*b[1]+a[1]*b[0];
}

__device__ __forceinline__ void complexAdd(const Number_t a[2], const Number_t b[2], Number_t out[2]){
	out[0] += a[0]+b[0];
	out[1] += a[1]+b[1];
}

__device__ __forceinline__ Number_t normalize_angle_2pi(Number_t ang){
	if(isnan(ang)) return 0;
	const Number_t PI = acos(-1.0f);
	ang = fmod(ang, 2*PI);
	if(ang < 0) ang += 2*PI;
	return ang;
}

__device__ __forceinline__ Number_t normalize_angle_pi(Number_t ang){
	const Number_t PI = acos(-1.0f);
	ang = fmod(ang+PI, 2*PI);
	if(ang < 0) ang += 2*PI;
	return ang-PI;
}

//K is the weight for a
//Assumes that a and b are in [-pi, pi]
__device__ __forceinline__ Number_t interpolate_angles(Number_t a, Number_t b, Number_t k){
	const Number_t PI = acos(-1.0f);
	Number_t diff = a-b;
	if(diff < 0){
		Number_t t = a;
		a = b;
		b = t;
		k = 1-k;
		diff = -diff;
	}
	diff = (diff > PI ? diff - 2*PI
					  : diff);

	return normalize_angle_pi(b + diff*k);
}