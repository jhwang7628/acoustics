#include "CUDA_PAT_WaveSolver.h"
#include <boost/bind.hpp>

#include <iostream>

Number_t gaussian_3d(const Number_t x, const Number_t y, const Number_t z, const Number_t b_x, const Number_t b_y, const Number_t b_z){
	Number_t stddev = 0.001;
	Number_t term = sqrt((x-b_x)*(x-b_x) + (y-b_y)*(y-b_y) + (z-b_z)*(z-b_z))/stddev;
	
	// Number_t term = x-mean;
	return exp(-term*term);
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

bool sphere(const Number_t x, const Number_t y, const Number_t z){
	Number_t dist = sqrt(x*x + y*y + z*z);
	if(dist < 0.01){
		return true;
	} else{
		return false;
	}
}

Number_t sphereNormal(const Number_t x, const Number_t y, const Number_t z, int dim){
	Number_t dist = sqrt(x*x + y*y + z*z);
	if(dim == 0){
		return x/dist;
	} else if(dim == 1){
		return y/dist;
	} else{
		return z/dist;
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
		Vector3d grad = field.gradient(v);
		grad.normalize();
		Number_t t = 0;
		int i, j, k;
		field.voxelIndices(vec, &i, &j, &k);
		ClosestPointField::Feature f;
		field.closestFeature(i, j, k, f);
		res += modedata.mode(mode)[3*f.index1 + dim]*f.alpha; t += f.alpha;
		res += modedata.mode(mode)[3*f.index2 + dim]*f.beta; t += f.beta;

		res += modedata.mode(mode)[3*f.index3 + dim]*f.gamma; t += f.gamma;
		if(t == 0) t = 1;
		res = res;
	}
	return res;
}

struct gambi_param{
	int nx;
	int ny;
	int nz;
	double dx;
	double dy;
	double dz;
	double xmin;
	double ymin;
	double zmin;
};

Number_t modeWithGambi(const Number_t x, const Number_t y, const Number_t z, int dim, const ClosestPointField & field, const std::vector<double> & dadata, const std::vector<int> & dacont, const gambi_param & params
						){
	int nx = params.nx;
	int ny = params.ny;
	int nz = params.nz;
	double dx = params.dx;
	double dy = params.dy;
	double dz = params.dz;
	double xmin = params.xmin;
	double ymin = params.ymin;
	double zmin = params.zmin;

	Vector3d v((REAL) x, (REAL) y, (REAL) z);
	Vec3d vec((REAL) x, (REAL) y, (REAL) z);
	Number_t res = 0;
	if(field.insideBox(vec)){
		Vector3d grad = field.gradient(v);
		grad.normalize();
		
		int xx = (int) floor((x - dx/2 - xmin)/dx);
		int yy = (int) floor((y - dy/2 - ymin)/dy);
		int zz = (int) floor((z - dz/2 - zmin)/dz);
		int pos = xx + nx*(yy + ny*zz);

		double prod = 0;
		int a = dacont[pos];
		if(a != 0){
			double xd = dadata[3*pos + 0]/a;
			double yd = dadata[3*pos + 1]/a;
			double zd = dadata[3*pos + 2]/a;
			prod = xd*grad[0] + yd*grad[1] + zd*grad[2];
		}

		res = grad[dim]*prod;
	}
	return res;
}

CUDA_PAT_WaveSolver::CUDA_PAT_WaveSolver(REAL timeStep,
										 const Vector3d & sound_source,
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
	this->_frequency = frequency/(2*acos(-1));
	_multipole_radius = multipole_radius;
	printf("<%lf, %lf, %lf>\n", xmax-xmin, ymax-ymin, zmax-zmin);
	printf("frequency: %f\n", frequency/(2*acos(-1)));
	
	_multipoleModeData._mode = -1;
	_multipoleModeData._frequency = frequency/(2*acos(-1));
	_multipoleModeData._coefficients.resize(2*(num_multipole_coef+1)*(num_multipole_coef+1));
	
	printf("%f <? %f\n", wave_speed*timeStep, cellSize);


	Number_t * posi = (Number_t *) malloc(3*listeningPositions->size()*sizeof(Number_t));

	for(int i = 0; i < listeningPositions->size(); i++){

		posi[3*i] = (Number_t) ((*listeningPositions)[i][0]);
		posi[3*i+1] = (Number_t) ((*listeningPositions)[i][1]);
		posi[3*i+2] = (Number_t) ((*listeningPositions)[i][2]);
	}

	// Wave_BoundaryEvaluator3D boundary = boost::bind(sphere, _1, _2, _3);
	// Wave_GradientEvaluator3D gradient = boost::bind(sphereNormal, _1, _2, _3, _4);
	Wave_InitialCondition3D initial = boost::bind(gaussian_3d, _1, _2, _3, sound_source[0], sound_source[1], sound_source[2]);
	// Wave_InitialCondition3D initial = boost::bind(gaussian_3d, _1, _2, _3, -0.0763, 0.003, -0.002);
	
	Wave_BoundaryEvaluator3D boundary = boost::bind(testWithDistanceField, _1, _2, _3, distanceTolerance, boost::ref(distanceField));
	Wave_GradientEvaluator3D gradient = boost::bind(gradientWithDistanceField, _1, _2, _3, _4, boost::ref(distanceField));
	

	this->wave = wave_sim_init(xmin, ymin, zmin,
							   xmax, ymax, zmax,
							   (Number_t) wave_speed, (Number_t) timeStep,
							   (Number_t) cellSize,
							   listeningPositions->size(),
							   posi,
							   initial,
							   boundary,
							   xcenter, ycenter, zcenter,
							   gradient,
							   (Number_t) pmlWidth*cellSize,
							   (Number_t) pmlStrength,
							   0,
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

#include "io/TetMeshReader.hpp"
#include <geometry/FixVtxTetMesh.hpp>
typedef FixVtxTetMesh<double>     TMesh;

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
	Number_t frequency = sqrt(modedata.omegaSquared(mode)/density)/(2*acos(-1));
	_multipole_radius = multipole_radius;

	_multipoleModeData._mode = mode;
	_multipoleModeData._frequency = frequency;
	_multipoleModeData._coefficients.resize(2*(num_multipole_coef+1)*(num_multipole_coef+1));

	this->_frequency = frequency;

	cellSize = min(cellSize, wave_speed/(8*frequency));
	printf("frequency: %f\n", frequency);
	printf("cellSize: %f\n", cellSize);
	timeStep = max(timeStep, cellSize/(2*wave_speed));
	printf("timeStep: %f\n", timeStep);
	printf("%f <? %f\n", wave_speed*timeStep, cellSize);

	Number_t * posi = (Number_t *) malloc(3*listeningPositions->size()*sizeof(Number_t));

	for(int i = 0; i < listeningPositions->size(); i++){

		posi[3*i] = (Number_t) ((*listeningPositions)[i][0]);
		posi[3*i+1] = (Number_t) ((*listeningPositions)[i][1]);
		posi[3*i+2] = (Number_t) ((*listeningPositions)[i][2]);
	}

	Wave_BoundaryEvaluator3D boundary = boost::bind(testWithDistanceField, _1, _2, _3, distanceTolerance, boost::ref(distanceField));
	
	/**
	 * HACK VIOLENTO
	 */
	const char * tmeshpath= "/home/rf356/Desktop/Code/acoustics/work/objects/ceramic_plate/ceramic_plate.tet";

	TMesh * tmesh = new TMesh();
	FV_TetMeshLoader_Double::load_mesh(tmeshpath, *tmesh);
	tmesh->init();
	tmesh->update_surface();

	int nx = ceil((xmax-xmin)/cellSize);
	int ny = ceil((ymax-ymin)/cellSize);
	int nz = ceil((zmax-zmin)/cellSize);
	double dx = (xmax-xmin)/nx;
	double dy = (ymax-ymin)/ny;
	double dz = (zmax-zmin)/nz;
	std::vector<double> modes(3*nx*ny*nz, 0.0);
	std::vector<int> cont(nx*ny*nz, 0.0);

	const std::vector<Point3d>& vs = tmesh->vertices();
	const std::vector<double> &damode = modedata.mode(mode);

	for(int i = 0; i < damode.size()/3; i++){
		int pos;
		int xx = (int) floor((vs[i].x - dx/2 - xmin)/dx);
		int yy = (int) floor((vs[i].y - dy/2 - ymin)/dy);
		int zz = (int) floor((vs[i].z - dz/2 - zmin)/dz);
		pos = xx + nx*(yy + ny*zz);
		modes[3*pos + 0] += damode[3*i + 0];
		modes[3*pos + 1] += damode[3*i + 1];
		modes[3*pos + 2] += damode[3*i + 2];
		cont[pos]++;
	}

	delete tmesh;

	gambi_param params;
	params.nx = nx;
	params.ny = ny;
	params.nz = nz;
	params.dx = dx;
	params.dy = dy;
	params.dz = dz;
	params.xmin = xmin;
	params.ymin = ymin;
	params.zmin = zmin;


	Wave_GradientEvaluator3D gradient = boost::bind(modeWithGambi, _1, _2, _3, _4, boost::ref(distanceField), boost::ref(modes), boost::ref(cont), boost::ref(params));
	// (const ClosestPointField & field, std::vector<double> & dadata, std::vector<int> & dacont,
	// 					int nx, int ny, int nz, double dx, double dy, double dz, double xmin, double ymin, double zmin)

	// Wave_GradientEvaluator3D gradient = boost::bind(modeWithClosestPointField, _1, _2, _3, _4, boost::ref(distanceField), mode, boost::ref(modedata));
	/*
	 * FIM DO HACK VIOLENTO	
	 */


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
							   frequency*2*acos(-1),
							   num_multipole_coef,
							   multipole_radius
							   );

	free(posi);

	// int nx, ny, nz;
	wave_sim_get_divisions(this->wave, &nx, &ny, &nz);
	printf("BOX: <%f - %f, %f - %f, %f - %f>\n", xmin, xmax, ymin, ymax, zmin, zmax);
	printf("BALL: (%f - %f, %f - %f, %f - %f)\n", xcenter-multipole_radius, xcenter+multipole_radius,
											ycenter-multipole_radius, ycenter+multipole_radius,
											zcenter-multipole_radius, zcenter+multipole_radius);
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
			this->_waveOutput[i][0].clear();
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
	if(_listeningPositions && (_step %     _substeps) == 0){
		Number_t * pressure = wave_listen(this->wave);
		for(int i = 0; i < _listeningPositions->size(); i++){
			_waveOutput[i][0].push_back((REAL) pressure[i]);
		}
	}


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
	Number_t x = wave_sim_get_x(wave, index[0]);
	Number_t y = wave_sim_get_y(wave, index[1]);
	Number_t z = wave_sim_get_z(wave, index[2]);
	Number_t r = sqrt(x*x + y*y + z*z);
	int pos = 4*(index[0] + nx*(index[1] + ny*index[2]));
	pressure[0] = (REAL) cache_pressure[pos];
	pos = index[0] + nx*(index[1] + ny*index[2]);
	pressure[1] = (REAL) cache_amplitude[pos];
	// Number_t k = 2*acos(-1)*4000/343.0;
	// pressure[2] = (0.01*0.01*1)/(r*sqrt(1+k*k*0.01*0.01))*cos(k*(r-0.01));
	pressure[2] = (REAL) cache_phase[pos];
	pressure[3] = (REAL) (cache_amplitude[pos]*cos(cache_phase[pos]));
}

void CUDA_PAT_WaveSolver::vertexData(int x, int y, int z, REAL * pressure, REAL * amplitude, REAL * phase, bool * bulk){
	VECTOR t(4);
	vertexPressure(Tuple3i(x, y, z), t);
	*pressure = t[0];
	*amplitude = t[1];
	*phase = t[2];
	*bulk = checkIsBulk(this->wave, x, y, z);
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

MultipoleData::MultipoleModeData & CUDA_PAT_WaveSolver::computeMultipoleCoefficients(){
	Number_t * coef = wave_compute_multipole(this->wave, _multipole_radius);
	for(int i = 0; i < _multipoleModeData._coefficients.size(); i++){
		_multipoleModeData._coefficients[i] = REAL(coef[i]);
	}
	//wave_test_multipole(this->wave);
	return _multipoleModeData;
}

REAL CUDA_PAT_WaveSolver::computeSound(const Vector3d & v, REAL t,  REAL * amplitude, REAL * freq,  REAL * phase){
	int i, j, k;
	wave_estimate_ijk(this->wave, v.x, v.y, v.z, &i, &j, &k);
	VECTOR p;
	Tuple3i tuple(i, j, k);
	this->vertexPressure(tuple, p);
	*amplitude = p[1];
	*freq = this->_frequency;
	*phase = p[2];
	return p[1]*cos(this->_frequency*t + p[2]);
}

REAL CUDA_PAT_WaveSolver::estimateSound(const Vector3d & v, REAL t,  REAL * amplitude, REAL * freq,  REAL * phase){
	Number_t aamplitude;
	Number_t aphase;
	wave_estimate_with_multipole(this->wave, v.x, v.y, v.z, &aamplitude, &aphase);
	*amplitude = aamplitude;
	*freq = this->_frequency;
	*phase = aphase;
	return aamplitude*cos(this->_frequency*t + aphase);
}

#include <multipole/MultipoleUtil.h>

void CUDA_PAT_WaveSolver::saveMultipoleCoefficients(const std::string & filename){
	computeMultipoleCoefficients();
	Multipole::saveToFile(filename, _multipoleModeData);
}

void CUDA_PAT_WaveSolver::gradientAt(int i, int j, int k, REAL * x, REAL * y, REAL * z){
	Number_t xx, yy, zz;
	wave_gradientAt(this->wave, i, j, k, &xx, &yy, &zz);
	*x = xx; *y = yy; *z = zz;
}