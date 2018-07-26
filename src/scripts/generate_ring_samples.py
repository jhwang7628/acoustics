#!/usr/bin/env python
import numpy as np

r = 0.135
offset = 0.0
center = [-0.2667, 0.6317, -0.1299]
N = 72
dim = 2

t = 0.0
dt = 2.0*np.pi/float(N)
for ii in range(N):
    vec = [0.0, 0.0, 0.0]
    vec[dim] = offset
    vec[(dim+1)%3] = r*np.cos(t)
    vec[(dim+2)%3] = r*np.sin(t)
    print '<listening_point x="%f" y="%f" z="%f"/>' %(vec[0], vec[1], vec[2])
    t += dt

