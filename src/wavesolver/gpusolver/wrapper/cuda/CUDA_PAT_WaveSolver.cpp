#include "CUDA_PAT_WaveSolver.h"
#include <boost/bind.hpp>

#include <iostream>

Number_t gaussian_3d(const Number_t x, const Number_t y, const Number_t z){
	Number_t stddev = 0.01*0.0568973;
	Number_t mean = -0.02275892;
	Number_t var2 = stddev*stddev*2;
	Number_t term = sqrt((x-mean)*(x-mean) + (y-mean)*(y-mean) + (z-0)*(z-0));
	// Number_t term = x-mean;
	return stddev*exp(-term*term/var2)/sqrt(acos(-1)*var2);
}

Number_t zeros(const Number_t x, const Number_t y, const Number_t z){
	return 0;
}

bool testWithDistanceField(const Number_t x, const Number_t y, const Number_t z, const REAL tolerance, const DistanceField & distanceField){
	Vector3d v((REAL) x, (REAL) y, (REAL) z);
	if(distanceField.distance(v) <= tolerance){
		return true;
	} else{
		return false;
	}
}

Number_t gradientWithDistanceField(const Number_t x, const Number_t y, const Number_t z, int dim, const DistanceField & distanceField){
	Vector3d v((REAL) x, (REAL) y, (REAL) z);
	Vector3d grad = distanceField.gradient(v);
	grad.normalize();
	return (Number_t) grad[dim];
}

Number_t modeWithClosestPointField(const Number_t x, const Number_t y, const Number_t z, int dim, const ClosestPointField & field, int mode, const ModeData & modedata){
	Vector3d v((REAL) x, (REAL) y, (REAL) z);
	Vec3d vec((REAL) x, (REAL) y, (REAL) z);
	Number_t res = 0;
	if(field.insideBox(vec)){
		Number_t t = 0;
		int i, j, k;
		field.voxelIndices(vec, &i, &j, &k);
		ClosestPointField::Feature f;
		field.closestFeature(i, j, k, f);
		res += modedata.mode(mode)[3*f.index1 + dim]*f.alpha; t += f.alpha;
		res += modedata.mode(mode)[3*f.index2 + dim]*f.beta; t += f.beta;

		res += modedata.mode(mode)[3*f.index3 + dim]*f.gamma; t += f.gamma;
		if(t == 0) t = 1;
		res = res/t;
	}
	return res;
}

CUDA_PAT_WaveSolver::CUDA_PAT_WaveSolver(REAL timeStep,
			        						   const BoundingBox & bbox, REAL cellSize,
			        						   const TriMesh & mesh,
			        						   const Vector3d & centerOfMass,
			        						   const DistanceField & distanceField,
			        						   REAL distanceTolerance,
			        						   const Vector3Array * listeningPositions,
			        						   WriteCallback * callback,
			        						   int substeps,
			        						   REAL endTime,
			        						   REAL frequency,
			        						   int num_multipole_coef,
			        						   REAL multipole_radius,
			        						   REAL pmlWidth,
			        						   REAL pmlStrength,
			        						   REAL wave_speed
			        						   ):_meshes(){

	_step = 0;
	cache_pressure = NULL;
    cache_amplitude = NULL;
    cache_phase = NULL;
	Number_t xmin = (Number_t) bbox.xmin();
	Number_t xmax = (Number_t) bbox.xmax();
	Number_t ymin = (Number_t) bbox.ymin();
	Number_t ymax = (Number_t) bbox.ymax();
	Number_t zmin = (Number_t) bbox.zmin();
	Number_t zmax = (Number_t) bbox.zmax();
	Number_t xcenter = (Number_t) centerOfMass[0];
	Number_t ycenter = (Number_t) centerOfMass[1];
	Number_t zcenter = (Number_t) centerOfMass[2];

	Number_t * posi = (Number_t *) malloc(3*listeningPositions->size()*sizeof(Number_t));

	for(int i = 0; i < listeningPositions->size(); i++){

		posi[3*i] = (Number_t) ((*listeningPositions)[i][0]);
		posi[3*i+1] = (Number_t) ((*listeningPositions)[i][1]);
		posi[3*i+2] = (Number_t) ((*listeningPositions)[i][2]);
	}

	Wave_BoundaryEvaluator3D boundary = boost::bind(testWithDistanceField, _1, _2, _3, distanceTolerance, boost::ref(distanceField));
	Wave_GradientEvaluator3D gradient = boost::bind(gradientWithDistanceField, _1, _2, _3, _4, boost::ref(distanceField));

	this->wave = wave_sim_init(xmin, ymin, zmin,
							   xmax, ymax, zmax,
							   (Number_t) wave_speed, (Number_t) timeStep,
							   (Number_t) cellSize,
							   listeningPositions->size(),
							   posi,
							   zeros,
							   boundary,
							   xcenter, ycenter, zcenter,
							   gradient,
							   (Number_t) pmlWidth*cellSize,
							   (Number_t) pmlStrength,
							   frequency,
							   num_multipole_coef,
							   multipole_radius
							   );

	free(posi);

	int nx, ny, nz;
	wave_sim_get_divisions(this->wave, &nx, &ny, &nz);

	this->_fieldDivisions = Vector3i(nx, ny, nz);
	this->_listeningPositions = listeningPositions;
	this->_meshes.push_back(&mesh);
	this->_callback = callback;
	this->_endTime = endTime;
	this->_substeps = substeps;
	if(listeningPositions){
		this->_waveOutput.resize(listeningPositions->size());
		for(int i = 0; i < listeningPositions->size(); i++){
			this->_waveOutput[i].resize(1);
		}
	}
}

CUDA_PAT_WaveSolver::CUDA_PAT_WaveSolver(REAL timeStep,
				                          const BoundingBox & bbox, REAL cellSize,
				                          const TriMesh & mesh,
				                          const Vector3d & centerOfMass,
				                          const ClosestPointField & distanceField,
				                          REAL distanceTolerance,
				                          int mode,
				                          const ModeData & modedata,
				                          REAL density,
				                          const Vector3Array * listeningPositions,
				                          WriteCallback * callback,
				                          int substeps,
				                          REAL endTime,
				                          int num_multipole_coef,
				                          REAL multipole_radius,
				                          REAL pmlWidth,
				                          REAL pmlStrength,
				                          REAL wave_speed
				                          ):_meshes(){

	_step = 0;
	cache_pressure = NULL;
    cache_amplitude = NULL;
    cache_phase = NULL;
	Number_t xmin = (Number_t) bbox.xmin();
	Number_t xmax = (Number_t) bbox.xmax();
	Number_t ymin = (Number_t) bbox.ymin();
	Number_t ymax = (Number_t) bbox.ymax();
	Number_t zmin = (Number_t) bbox.zmin();
	Number_t zmax = (Number_t) bbox.zmax();
	Number_t xcenter = (Number_t) centerOfMass[0];
	Number_t ycenter = (Number_t) centerOfMass[1];
	Number_t zcenter = (Number_t) centerOfMass[2];
	Number_t frequency = sqrt(modedata.omegaSquared(mode)/density);

	printf("frequency: %f\n", frequency/(2*acos(-1)));

	Number_t * posi = (Number_t *) malloc(3*listeningPositions->size()*sizeof(Number_t));

	for(int i = 0; i < listeningPositions->size(); i++){

		posi[3*i] = (Number_t) ((*listeningPositions)[i][0]);
		posi[3*i+1] = (Number_t) ((*listeningPositions)[i][1]);
		posi[3*i+2] = (Number_t) ((*listeningPositions)[i][2]);
	}

	Wave_BoundaryEvaluator3D boundary = boost::bind(testWithDistanceField, _1, _2, _3, distanceTolerance, boost::ref(distanceField));
	Wave_GradientEvaluator3D gradient = boost::bind(modeWithClosestPointField, _1, _2, _3, _4, boost::ref(distanceField), mode, boost::ref(modedata));

	this->wave = wave_sim_init(xmin, ymin, zmin,
							   xmax, ymax, zmax,
							   (Number_t) wave_speed, (Number_t) timeStep,
							   (Number_t) cellSize,
							   listeningPositions->size(),
							   posi,
							   zeros,
							   boundary,
							   xcenter, ycenter, zcenter,
							   gradient,
							   (Number_t) pmlWidth*cellSize,
							   (Number_t) pmlStrength,
							   frequency,
							   num_multipole_coef,
							   multipole_radius
							   );

	free(posi);

	int nx, ny, nz;
	wave_sim_get_divisions(this->wave, &nx, &ny, &nz);

	this->_fieldDivisions = Vector3i(nx, ny, nz);
	this->_listeningPositions = listeningPositions;
	this->_meshes.push_back(&mesh);
	this->_callback = callback;
	this->_endTime = endTime;
	this->_substeps = substeps;
	if(listeningPositions){
		this->_waveOutput.resize(listeningPositions->size());
		for(int i = 0; i < listeningPositions->size(); i++){
			this->_waveOutput[i].resize(1);
		}
	}
}

CUDA_PAT_WaveSolver::~CUDA_PAT_WaveSolver(){
	wave_sim_free(this->wave);
}

const Tuple3i & CUDA_PAT_WaveSolver::fieldDivisions() const{
	return this->_fieldDivisions;
}

bool CUDA_PAT_WaveSolver::stepSystem(const BoundaryEvaluator &bcEvaluator){
	wave_sim_step(this->wave);
	REAL time = (REAL) wave_sim_get_current_time(this->wave);

	
	//Save output
	// if(_listeningPositions && (_step % _substeps) == 0){
	// 	for(int field = 0; field < 6; field++){
	// 		Number_t * pressure = wave_listen(this->wave, field);
	// 		for(int i = 0; i < _listeningPositions->size(); i++){
	// 			_waveOutput[i][field].push_back((REAL) pressure[i]);
	// 		}
	// 	}
	// }


	_step++;
	if(this->_endTime > 0 && time > this->_endTime){
		return false;
	}
	return true;
}

//TODO
void CUDA_PAT_WaveSolver::writeWaveOutput() const{
	if(!this->_callback){
		return;
	}
	(*(this->_callback))(this->_waveOutput);
}

Vector3d CUDA_PAT_WaveSolver::fieldPosition(const Tuple3i & index) const{
	REAL x = (REAL) wave_sim_get_x(this->wave, index[0]);
	REAL y = (REAL) wave_sim_get_y(this->wave, index[1]);
	REAL z = (REAL) wave_sim_get_z(this->wave, index[2]);

	return Vector3d(x, y, z);
}

Vector3d CUDA_PAT_WaveSolver::fieldPosition(int index) const{
	int nx, ny, nz;
	wave_sim_get_divisions(this->wave, &nx, &ny, &nz);
	int id = index;
	int i = index % (nx*ny);
	id = (id-i)/nx;
	int j = id % ny;
	int k = (id - j)/ny;

	REAL x = (REAL) wave_sim_get_x(this->wave, i);
	REAL y = (REAL) wave_sim_get_y(this->wave, j);
	REAL z = (REAL) wave_sim_get_z(this->wave, k);

	return Vector3d(x, y, z);
}

const Vector3Array * CUDA_PAT_WaveSolver::listeningPositions() const{
	return this->_listeningPositions;
}

//EXTREMELY SLOW
//USE ONLY TO DEBUG
void CUDA_PAT_WaveSolver::vertexPressure(const Tuple3i & index,
										 VECTOR & pressure){

	if(cache_pressure == NULL || _step % _substeps == 0){
		cache_pressure = wave_sim_get_u(this->wave);
	}

	if(cache_amplitude == NULL || _step % _substeps == 0){
		cache_amplitude = wave_sim_get_amplitudes(this->wave);
	}

	if(cache_phase == NULL || _step % _substeps == 0){
		cache_phase = wave_sim_get_phases(this->wave);
	}

	if(pressure.size() != 4){
		pressure.resizeAndWipe(4);
	}

	int nx, ny, nz;
	wave_sim_get_divisions(this->wave, &nx, &ny, &nz);
	int pos = 4*(index[0] + nx*(index[1] + ny*index[2]));
	pressure[0] = (REAL) cache_pressure[pos];
	pos = index[0] + nx*(index[1] + ny*index[2]);
	pressure[1] = (REAL) cache_amplitude[pos];
	pressure[2] = (REAL) cache_phase[pos];
	pressure[3] = (REAL) (cache_amplitude[pos]*cos(cache_phase[pos]));
}

int CUDA_PAT_WaveSolver::numCells() const{
	int nx, ny, nz;
	wave_sim_get_divisions(this->wave, &nx, &ny, &nz);
	return nx*ny*nz;
}

REAL CUDA_PAT_WaveSolver::currentSimTime() const{
	REAL time = (REAL) wave_sim_get_current_time(this->wave);
	return time;
}

REAL CUDA_PAT_WaveSolver::fieldDiameter() const{
	Number_t xmin, xmax, ymin, ymax, zmin, zmax;
	wave_sim_get_bounds(this->wave, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);
	Number_t ret = xmax-xmin;
	if(ymax-ymin > ret) ret = ymax-ymin;
	if(zmax-zmin > ret) ret = zmax-zmin;

	return (REAL)ret;
}

Vector3d CUDA_PAT_WaveSolver::sceneCenter() const{
	Number_t xmin, xmax, ymin, ymax, zmin, zmax;
	wave_sim_get_bounds(this->wave, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);
	return Vector3d((REAL)((xmax+xmin)/2), (REAL)((ymax+ymin)/2), (REAL)((zmax+zmin)/2));
}

void CUDA_PAT_WaveSolver::computeMultipoleCoefficients(){
	wave_test_multipole(this->wave);
}