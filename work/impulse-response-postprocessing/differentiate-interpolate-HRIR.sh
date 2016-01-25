#/bin/bash

BIN=/home/jui-hsien/code/acoustics/build/bin

# ${BIN}/differentiate-interpolate-HRIR /home/jui-hsien/code/acoustics/work/impulse-response/test_freespace/test_set head_left_pressure_ 100 100 100 -1.701 -1.701 -1.701 1.701 1.701 1.701 /home/jui-hsien/code/turbsound_postfluid/FOAM_templates/smooth_BL2_1m_refineconcave/cellCentroid.dat HRTFcatedR.dat
# ${BIN}/differentiate-interpolate-HRIR /home/jui-hsien/code/acoustics/work/impulse-response/test_freespace/test_set head_left_pressure_ 100 100 100 -1.701 -1.701 -1.701 1.701 1.701 1.701 /home/jui-hsien/code/turbsound_postfluid/FOAM_templates/smooth_BL2_1m_refineconcave/cellCentroid.dat HRTFcatedR.dat
${BIN}/differentiate-interpolate-HRIR /home/jui-hsien/code/acoustics/work/impulse-response/data_Gaussian_space_1cell_time_1over1E5_res250_simulation/test_set head_left_pressure_ 250 250 250  -0.2178337591671500 -0.2017749603171500 -0.2192579190471500 0.2204221281671500 0.2364809270171500 0.2189979682871500 /home/jui-hsien/code/turbsound_postfluid/FOAM_templates/smooth_BL2_1m_refineconcave/cellCentroid.dat HRTFcatedL.dat
