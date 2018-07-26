#!/usr/bin/env python
import numpy as np

h = 0.002
N = 150 - 6*2
center = np.array([-0.2667, 0.6317, -0.1299])
corner = center - N/2.*h

h = 0.01
N = 28

pts = []
for kk in range(N):
    for jj in range(N):
        for ii in range(N):
            if (ii==0 or ii==N-1 or
                jj==0 or jj==N-1 or
                kk==0 or kk==N-1):
                pos = corner.copy()
                pos[0] += (0.5+ii)*h
                pos[1] += (0.5+jj)*h
                pos[2] += (0.5+kk)*h
                pts.append(pos)
                print '<listening_point x="%f" y="%f" z="%f"/>' %(pos[0], pos[1], pos[2])
# print pts
