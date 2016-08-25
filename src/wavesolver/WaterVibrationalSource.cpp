#include <wavesolver/WaterVibrationalSource.h> 
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <sndgen/WavReader.hpp>
#include <utils/STL_Wrapper.h>
#include <fstream> 

//##############################################################################
//##############################################################################
WaterVibrationalSource::
WaterVibrationalSource(RigidObjectPtr owner, const std::string &wavFile)
    : VibrationalSource(owner), _surfaceMesh(owner->GetMeshPtr())
{
    Initialize(wavFile); 
    // FIXME debug
    TestSpline(); 
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    return 0.0; 
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSource::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented"); 
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSource::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented"); 
}

//##############################################################################
//##############################################################################
void WaterVibrationalSource::
Initialize(const std::string &wavFile)
{
    std::cout << "Initialize WaterVibrationalSource with file: " << wavFile << std::endl;
    ReadOscillatorFromWav(wavFile); 
    ComputeVelocityAndAcceleration();
    PrecomputeInterpolation(); 
}

//##############################################################################
//##############################################################################
void WaterVibrationalSource::
ReadOscillatorFromWav(const std::string &wavFile)
{
    std::cout << " Read oscillator displacement from wav file\n"; 

    // prepare and read from hydrophone channel (default to be right channel: 1)
    WavReader<REAL> reader; 
    reader.Open(wavFile); 
    reader.ReadChannel(_oscillatorDisplacement, 1);
    STL_Wrapper::PrintVectorContent(std::cout, _oscillatorDisplacement, 10);
    _sampleRate = reader.SampleRate(); 
    reader.Close(); 
}

//##############################################################################
// This function computes velocity and acceleration at given points. 
//
// Interior points always have 2nd accuracy. At the boundary, velocity has 2nd
// accuracy and acceleration has 1st accuracy. Reference for the FD
// coefficients: 
//  https://en.wikipedia.org/wiki/Finite_difference_coefficient
//##############################################################################
void WaterVibrationalSource::
ComputeVelocityAndAcceleration()
{
    std::cout << " Compute discrete velocity and acceleration signal from oscillator displacement\n"; 

    // make sure oscillator displacement is read and initialized
    const int N_frames = _oscillatorDisplacement.size(); 
    assert(N_frames > 0); 
    _oscillatorVelocity.resize(N_frames); 
    _oscillatorAcceleration.resize(N_frames); 

    const REAL half_over_h = 0.5 / _sampleRate; 
    const REAL one_over_h  = 1.0 / _sampleRate; 
    const REAL one_over_h2 = 1.0 / pow(_sampleRate, 2); 
    // interior point finite difference
    for (int f_idx=1; f_idx<N_frames-1; ++f_idx)
    {
        const REAL &v_f = _oscillatorDisplacement.at(f_idx+1); 
        const REAL &v_b = _oscillatorDisplacement.at(f_idx-1); 
        const REAL &v_c = _oscillatorDisplacement.at(f_idx); 
        _oscillatorVelocity.at(f_idx) = (v_f - v_b) * half_over_h; 
        _oscillatorAcceleration.at(f_idx) = (v_b + v_f - v_c) * one_over_h2; 
    }

    // first point
    {
        const int f_idx = 0;
        const REAL &v_0 = _oscillatorDisplacement.at(f_idx); 
        const REAL &v_1 = _oscillatorDisplacement.at(f_idx+1); 
        const REAL &v_2 = _oscillatorDisplacement.at(f_idx+2); 
        _oscillatorVelocity.at(f_idx) = (-1.5*v_0 + 2.0*v_1 - 1.5*v_2) * one_over_h;   
        _oscillatorAcceleration.at(f_idx) = (v_0 -2.0*v_1 + v_2) * one_over_h2;
    }

    // last point (note the sign)
    {
        const int f_idx = N_frames - 1;
        const REAL &v_0 = _oscillatorDisplacement.at(f_idx); 
        const REAL &v_1 = _oscillatorDisplacement.at(f_idx-1); 
        const REAL &v_2 = _oscillatorDisplacement.at(f_idx-2); 
        _oscillatorVelocity.at(f_idx) = (-1.5*v_0 + 2.0*v_1 - 1.5*v_2) * (-1.0*one_over_h);   
        _oscillatorAcceleration.at(f_idx) = (v_0 -2.0*v_1 + v_2) * (+1.0*one_over_h2);
    }

    STL_Wrapper::PrintVectorContent(std::cout, _oscillatorDisplacement, 10); 
    STL_Wrapper::PrintVectorContent(std::cout, _oscillatorVelocity, 10); 
    STL_Wrapper::PrintVectorContent(std::cout, _oscillatorAcceleration, 10); 
}

//##############################################################################
// This function initializes the resources for velocity and accleration
// interpolation.
//##############################################################################
void WaterVibrationalSource::
PrecomputeInterpolation()
{
    std::cout << "Precompute interpolator for velocity and acceleration signals.\n";

    assert(_oscillatorVelocity.size()>0 && _oscillatorAcceleration.size()>0); 
    const int N_frames = _oscillatorAcceleration.size(); 
    _oscillatorTime.resize(N_frames); 
    //FloatArray oscillatorTime(N_frames); 
    for (int t_idx=0; t_idx<N_frames; ++t_idx)
        _oscillatorTime.at(t_idx) = _startTime + _sampleRate*(REAL)t_idx; 
    _interpolatorAcceleration.init(N_frames, &_oscillatorTime[0], &_oscillatorAcceleration[0]); 
    std::cout << " Preocompute completed.\n"; 
}

//##############################################################################
//##############################################################################
void WaterVibrationalSource::
TestSpline()
{
    std::cout << "Testing acceleration interpolator C-Spline\n"; 

    std::ofstream ofRaw("rawSignal.txt"); 
    std::ofstream ofInterp("interpolatedSignal.txt"); 

    std::cout << " Interpolator Time range:\n"; 
    STL_Wrapper::PrintVectorContent(std::cout, _oscillatorTime, 10); 
    const int N_frames = _oscillatorAcceleration.size(); 
    for (int f_idx=0; f_idx<N_frames; ++f_idx)
        ofRaw << _oscillatorTime.at(f_idx) << " " << _oscillatorAcceleration.at(f_idx) << std::endl;
    ofRaw.close(); 


    const int factor = 10;
    const int N_interpFrames = N_frames * factor - (factor-1); // end at last frame of raw
    const REAL interpSampleRate = _sampleRate / (REAL)factor; 
    FloatArray interpTime; 
    for (int f_idx=0; f_idx<N_interpFrames; ++f_idx)
    {
        const REAL time = _startTime + (REAL)f_idx*interpSampleRate; 
        interpTime.push_back(time); 
        ofInterp << time << " " << _interpolatorAcceleration.eval(time) << std::endl;
    }
    ofInterp.close(); 
    std::cout << " Queried Time range:\n"; 
    STL_Wrapper::PrintVectorContent(std::cout, interpTime, 10);

    std::cout << "Results written to files 'rawSignal.txt' and 'interpolatedSignal.txt'" << std::endl;
}
