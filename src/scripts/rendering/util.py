#!/usr/bin/env python
def Resample(step_freq, target_freq, target_total_frames):
    assert target_freq < step_freq
    skip_frames = int(float(step_freq)/target_freq)
