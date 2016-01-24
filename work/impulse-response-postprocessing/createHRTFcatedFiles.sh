#/bin/bash

./InterpolateHRIR /home/jui-hsien/code/acoustics/work/impulse-response/data_Gaussian0p0025_res250_offset2stddev_simulation head_right_pressure_ 250 250 250 -0.2178337591671500 -0.2017749603171500 -0.2192579190471500 0.2204221281671500 0.2364809270171500 0.2189979682871500 /home/jui-hsien/code/turbsound_postfluid/FOAM_templates/smooth_BL2_1m_refineconcave/cellCentroid.dat HRTFcatedR.dat
./InterpolateHRIR /home/jui-hsien/code/acoustics/work/impulse-response/data_Gaussian0p0025_res250_offset2stddev_simulation head_left_pressure_ 250 250 250  -0.2178337591671500 -0.2017749603171500 -0.2192579190471500 0.2204221281671500 0.2364809270171500 0.2189979682871500 /home/jui-hsien/code/turbsound_postfluid/FOAM_templates/smooth_BL2_1m_refineconcave/cellCentroid.dat HRTFcatedL.dat
