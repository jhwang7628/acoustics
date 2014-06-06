#ifndef CUDA_PAN_SOLVER_H
#define CUDA_PAN_SOLVER_H

#include "../../cuda/cuda_PAN_wave_3d.h"
#include <distancefield/distanceField.h>

#include <geometry/BoundingBox.h>
#include <geometry/TriangleMesh.hpp>

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>
#include <linearalgebra/Vector3.hpp>

#include <TYPES.h>
#include "../../../WaveSolver.h"

#include <utils/Evaluator.h>

#include <vector>

class CUDA_PAN_WaveSolver : public Solver {
    public:
    	typedef boost::function<void (const vector<vector<FloatArray> >&w)> WriteCallback;

    public:
    	CUDA_PAN_WaveSolver(REAL timeStep,
        						  const BoundingBox & bbox, REAL cellSize,
        						  const TriMesh & mesh,
        						  const Vector3d & centerOfMass,
        						  const DistanceField & distanceField,
        						  REAL distanceTolerance,
        						  const Vector3Array * listeningPositions,
        						  WriteCallback * callback,
        						  int substeps,
        						  REAL endTime,
        						  REAL pulseTime,
        						  REAL pmlWidth=11.0,
        						  REAL pmlStrength=1000000.0,
        						  REAL wave_speed=343.0
        						  );

        virtual ~CUDA_PAN_WaveSolver();

        virtual int N() const{
        	return 6;
        }

        virtual const Tuple3i &fieldDivisions() const;
        virtual const std::vector<const TriMesh *> &meshes() const{
            return this->_meshes;
        }
        virtual bool stepSystem( const BoundaryEvaluator &bcEvaluator );
        virtual void writeWaveOutput() const;
        virtual Vector3d fieldPosition( const Tuple3i &index ) const;
        virtual Vector3d fieldPosition( int index ) const;
        virtual const Vector3Array *listeningPositions() const;
        virtual REAL fieldDiameter() const;
        virtual void vertexPressure( const Tuple3i &index, VECTOR &pressure );
        virtual int numCells() const;
        virtual REAL currentSimTime() const;
        virtual Vector3d sceneCenter() const;

    private:
    	Cuda_PAN_Wave_3d_t wave;

    	const Vector3Array * _listeningPositions;

    	std::vector<const TriMesh *> _meshes;
    	std::vector<std::vector<FloatArray> > _waveOutput;
    	WriteCallback * _callback;
    	REAL _endTime;
    	int _substeps;
        int _step;
        Vector3i _fieldDivisions;
        Number_t * cache;
};

#endif