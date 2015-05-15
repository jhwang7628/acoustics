#ifndef CUDA_PAT_SOLVER_H
#define CUDA_PAT_SOLVER_H

#include "../../cuda/cuda_PAT_wave_3d.h"
#include <distancefield/distanceField.h>
#include <distancefield/closestPointField.h>

#include <deformable/ModeData.h>

#include <geometry/BoundingBox.h>
#include <geometry/TriangleMesh.hpp>

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>
#include <linearalgebra/Vector3.hpp>

#include <TYPES.h>
#include "../../../WaveSolver.h"

#include <utils/Evaluator.h>

#include <vector>
#include <string>

#include <multipole/MultipoleUtil.h>

class CUDA_PAT_WaveSolver : public Solver {
		public:
			typedef boost::function<void (const vector<vector<FloatArray> >&w)> WriteCallback;

		public:

			//Needs refactor
			CUDA_PAT_WaveSolver(REAL timeStep, // this is used in precompute_impulse_response_text
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
								REAL pmlWidth=15.0,
								REAL pmlStrength=1000000.0,
								REAL wave_speed=343.0
								);

			CUDA_PAT_WaveSolver(REAL timeStep,
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
								REAL pmlWidth=15.0,
								REAL pmlStrength=1000000.0,
								REAL wave_speed=343.0
								);

			virtual ~CUDA_PAT_WaveSolver();

			virtual int N() const{
					//For debug only
					//0 - pressure
					//1 - Amplitude
					//2 - Phase
					//3 - Amplitude*cos(Phase)
				return 4;
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
			virtual MultipoleData::MultipoleModeData & computeMultipoleCoefficients();
			virtual REAL computeSound(const Vector3d & v, REAL t,  REAL * amplitude, REAL * freq,  REAL * phase);
			virtual REAL estimateSound(const Vector3d & v, REAL t,  REAL * amplitude, REAL * freq,  REAL * phase);
			virtual void saveMultipoleCoefficients(const std::string & filename);

			void vertexData(int x, int y, int z, REAL * pressure, REAL * amplitude, REAL * phase, bool * bulk);
			void gradientAt(int i, int j, int k, REAL * x, REAL * y, REAL * z);

		private:
			Cuda_PAT_Wave_3d_t wave;

			const Vector3Array * _listeningPositions;
			REAL _frequency;
			std::vector<const TriMesh *> _meshes;
			std::vector<std::vector<FloatArray> > _waveOutput;
			MultipoleData::MultipoleModeData _multipoleModeData;
			WriteCallback * _callback;
			REAL _endTime;
			REAL _multipole_radius;
			int _substeps;
			int _step;
			Vector3i _fieldDivisions;
			Number_t * cache_pressure;
			Number_t * cache_amplitude;
			Number_t * cache_phase;
};

#endif
