#ifndef GPU_WAVE_DEFINES_H
#define GPU_WAVE_DEFINES_H

#include <boost/function.hpp>

	typedef float Number_t;

	typedef boost::function<Number_t (const Number_t x, const Number_t y, const Number_t z)>
		Wave_InitialCondition3D;
	typedef boost::function<bool (const Number_t x, const Number_t y, const Number_t z)>
		Wave_BoundaryEvaluator3D;
	typedef boost::function<Number_t (const Number_t x, const Number_t y, const Number_t z, int dim)>
		Wave_GradientEvaluator3D;

#endif