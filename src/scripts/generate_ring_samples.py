#!/usr/bin/env python
import numpy as np

r = 1.0
offset = 0.716
N = 36

t = 0.0
dt = 2.0*np.pi/float(N)
for ii in range(N):
    x = r*np.cos(t)
    z = r*np.sin(t)
    y = offset
    print '<listening_point x="%f" y="%f" z="%f"/>' %(x,y,z)
    t += dt

