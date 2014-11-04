//////////////////////////////////////////////////////////////////////
// precompute_acceleration_pulse.cpp: Precomputes the result of a
//                                    body accelerating over a short
//                                    time scale where the
//                                    acceleration pulse is modelled
//                                    using a simple interpolation
//                                    function
//
//////////////////////////////////////////////////////////////////////

#include <QtGui>

#include <config.h>
#include <TYPES.h>

#include <distancefield/closestPointField.h>
#include <distancefield/FieldBuilder.h>

#include <geometry/RigidMesh.h>
#include <geometry/TriangleMesh.hpp>

#include <deformable/ModeData.h>

#include <linearalgebra/Vector3.hpp>

#include <math/InterpolationFunction.h>

#include <parser/Parser.h>

#include <transfer/PulseApproximation.h>

#include <ui/WaveViewer.h>

#include <utils/IO.h>
#include <utils/MathUtil.h>

#ifdef USE_CUDA
	#include <wavesolver/gpusolver/wrapper/cuda/CUDA_PAN_WaveSolver.h>
	#include <wavesolver/gpusolver/wrapper/cuda/CUDA_PAT_WaveSolver.h>
#else
	#include <wavesolver/PML_WaveSolver.h>
#endif

#include <wavesolver/WaveSolver.h>

#include <boost/bind.hpp>

#include <QGLViewer/qglviewer.h>

#include <GL/glut.h>

#include <iostream>
#include <string>

using namespace std;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
#include <boost/python.hpp>

using namespace boost::python;

class PAT_Wrapper{

public:
	Parser::AcousticTransferParms parms;
	Vector3Array             listeningPositions;
	TriangleMesh<REAL> * mesh;
	RigidMesh * rigidMesh;
	std::string meshFileName;
	CUDA_PAT_WaveSolver * solver;

	Vector3d centerOfMass;
	
	ClosestPointField * sdf;

	BoundingBox fieldBBox;

	//To Python
	ModeData modeData;
	REAL radius;
	REAL density;

	REAL cellSize;

	REAL timeStep;
	int sdfResolution;
	REAL gridScale;
	int gridResolution;
	int mode;
	int nbar;
	int substeps;
	REAL endTime;
	REAL scaleRadius;
	REAL wave_speed;

public:
	PAT_Wrapper(const std::string & fileName="default.xml"){
		rigidMesh = NULL;
		mesh = NULL;
		solver = NULL;
		sdf = NULL;
		wave_speed = 343.0;
		Parser * parser = Parser::buildParser( fileName );

		if ( !parser )
		{
			cerr << "ERROR: Could not build parser from " << fileName << endl;
		}

		mesh = parser->getMesh();

		if ( !mesh )
		{
			cerr << "ERROR: Could not build mesh" << endl;
		}

		meshFileName = parser->getMeshFileName();
		parms = parser->getAcousticTransferParms();

		rigidMesh = new RigidMesh( *mesh, parms._rigidPrefix, parms._rigidDensity );


		centerOfMass = rigidMesh->centerOfMass();

		//Data
		radius = mesh->boundingSphereRadius( rigidMesh->centerOfMass() );
		modeData.read(parms._modeDataFile.c_str());
		density = parms._rigidDensity;

		//Parameters
		timeStep = 1.0 / (REAL)( parms._timeStepFrequency);
		sdfResolution = parms._sdfResolution;
		gridScale = parms._gridScale;
		gridResolution = parms._gridResolution;
		mode = parms._mode;
		nbar = parms._nbar;
		substeps = parms._subSteps;
		endTime = -1;
		scaleRadius = parms._radiusMultipole;
	}

	~PAT_Wrapper(){
		if(rigidMesh){ delete rigidMesh;}
		if(mesh){ delete mesh;}
		if(solver){ delete solver;}
		if(sdf){ delete sdf;}
	}
	
public:
	void initSDF(){
		sdf = DistanceFieldBuilder::BuildSignedClosestPointField(meshFileName.c_str(),
															  sdfResolution,
															  parms._sdfFilePrefix.c_str() );
		fieldBBox = BoundingBox( sdf->bmin(), sdf->bmax() );

		// Scale this up to build a finite difference field
		fieldBBox *= gridScale;

		int cellDivisions = gridResolution;

		cellSize = min( fieldBBox.axislength( 0 ),
				min( fieldBBox.axislength( 1 ), fieldBBox.axislength( 2 ) ) );
		
		cellSize /= (REAL)cellDivisions;
	}

	void initSolver(){
		solver = new CUDA_PAT_WaveSolver(
				   timeStep,
				   fieldBBox, cellSize,
				   *mesh, centerOfMass,
				   *sdf,
				   0.0,
				   mode, //Mode
				   modeData,
				   density,
				   &listeningPositions, //listeningPositions
				   NULL,
				   substeps,
				   endTime,
				   nbar,
				   radius*scaleRadius,
				   50,
				   100000);
	}

	void setEndTime(double alpha){
		endTime =  alpha*timeStep + radius*scaleRadius/wave_speed;
	}

	void stepSolver(){
		solver->stepSystem(NULL);
	}

	void runSolver(){
		bool keep = true;
		while(keep){
			keep = solver->stepSystem(NULL);
		}
	}

	void saveToFile(const std::string & file){
		solver->saveMultipoleCoefficients(file);
	}

	tuple cellPosition(int x, int y, int z){
		Vector3d pos = solver->fieldPosition(Tuple3i(x, y, z));
		return make_tuple(pos.x, pos.y, pos.z);
	}

	tuple cellData(int x, int y, int z){
		double pressure, amplitude, phase;
		bool bulk;
		solver->vertexData(x, y, z, &pressure, &amplitude, &phase, &bulk);
		return make_tuple(pressure, std::complex<double>(amplitude*cos(phase), amplitude*sin(phase)), bulk);
	}

	tuple cellDivisions(){
		Tuple3i t = solver->fieldDivisions();
		return make_tuple(t[0], t[1], t[2]);
	}

	tuple gradientAt(int x, int y, int z){
		double xx, yy, zz;
		solver->gradientAt(x, y, z, &xx, &yy, &zz);
		return make_tuple(xx, yy, zz);
	}

};

BOOST_PYTHON_MODULE(_solver)
{
	class_<PAT_Wrapper>("PAT_Solver", init<std::string>())
		.def_readwrite("radius", &PAT_Wrapper::radius)
		.def_readwrite("density", &PAT_Wrapper::density)
		.def_readwrite("cellSize", &PAT_Wrapper::cellSize)
		.def_readwrite("timeStep", &PAT_Wrapper::timeStep)
		.def_readwrite("sdfResolution", &PAT_Wrapper::sdfResolution)
		.def_readwrite("gridScale", &PAT_Wrapper::gridScale)
		.def_readwrite("gridResolution", &PAT_Wrapper::gridResolution)
		.def_readwrite("mode", &PAT_Wrapper::mode)
		.def_readwrite("nbar", &PAT_Wrapper::nbar)
		.def_readwrite("substeps", &PAT_Wrapper::substeps)
		.def_readwrite("endTime", &PAT_Wrapper::endTime)
		.def_readwrite("scaleRadius", &PAT_Wrapper::scaleRadius)
		.def_readwrite("wave_speed", &PAT_Wrapper::wave_speed)
		.def_readwrite("modeData", &PAT_Wrapper::modeData)
		.def("initSDF", &PAT_Wrapper::initSDF)
		.def("initSolver", &PAT_Wrapper::initSolver)
		.def("setEndTime", &PAT_Wrapper::setEndTime)
		.def("stepSolver", &PAT_Wrapper::stepSolver)
		.def("runSolver", &PAT_Wrapper::runSolver)
		.def("saveToFile", &PAT_Wrapper::saveToFile)
		.def("cellPosition", &PAT_Wrapper::cellPosition)
		.def("cellData", &PAT_Wrapper::cellData)
		.def("cellDivisions", &PAT_Wrapper::cellDivisions)
		.def("gradientAt", &PAT_Wrapper::gradientAt);
}