#!/usr/bin/env python 
import math
import numpy as np 

def GenerateLine(x, y, z): 
    return '<listening_point x=\"%f\" y=\"%f\" z=\"%f\"/>' %(x, y, z)

def SampleCircle_Z(center, radius, N): 
    print 'Sampling on a z-circle'
    dt = 2.0*math.pi/float(N)
    for idx in range(N): 
        theta = dt*idx
        x = radius * math.cos(theta)
        y = radius * math.sin(theta)
        z = center
        print GenerateLine(x, y, z)

def SampleLine(start, end, N): 
    print 'Sampling on a line'
    L = np.linalg.norm(end - start)
    dl = L/float(N-1)
    normal = (end - start) / L
    for idx in range(N): 
        point = start + normal * (dl*idx)
        print GenerateLine(point[0], point[1], point[2])

if __name__ == '__main__': 
    # SampleCircle_Z(0.0, 0.15, 200)
    SampleLine(np.array([0.05, 0.0, 0.0]), np.array([0.2, 0.0, 0.0]), 20)
