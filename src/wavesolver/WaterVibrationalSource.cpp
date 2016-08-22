#include <wavesolver/WaterVibrationalSource.h> 
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <sndgen/WavReader.hpp>
#include <utils/STL_Wrapper.h>

//##############################################################################
//##############################################################################
WaterVibrationalSource::
WaterVibrationalSource(RigidObjectPtr owner, const std::string &wavFile)
    : VibrationalSource(owner), _surfaceMesh(owner->GetMeshPtr())
{
    Initialize(wavFile); 
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
    _oscillatorSampleRate = reader.SampleRate(); 
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

    const REAL half_over_h = 0.5 / _oscillatorSampleRate; 
    const REAL one_over_h  = 1.0 / _oscillatorSampleRate; 
    const REAL one_over_h2 = 1.0 / pow(_oscillatorSampleRate, 2); 
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
