#!/usr/bin/env python 
import math
import numpy as np 

def ObjectPlacementTemplate(x, y, z, ID):
    rigidObjectStr = """
            <rigid_object 
                working_directory="/home/jui-hsien/code/acoustics/work/two_sphere_collide"
                object_prefix="sphere_d10cm"  
                id="%u" 
                fieldresolution="400"
                scale="1.0" 
                initial_position_x="%.8f"
                initial_position_y="%.8f"
                initial_position_z="%.8f"
                material_id="0"
                /> """ %(ID, x, y, z)

    print rigidObjectStr
    return rigidObjectStr

def GenerateLine(x, y, z): 
    return '<listening_point x=\"%f\" y=\"%f\" z=\"%f\"/>' %(x, y, z)

def SampleCircle_Z(center, radius, N): 
    dt = 2.0*math.pi/float(N)
    for idx in range(N): 
        theta = dt*idx
        x = radius * math.cos(theta)
        y = radius * math.sin(theta)
        z = center
        print GenerateLine(x, y, z)

def SampleCircle_X(center, radius, N): 
    dt = 2.0*math.pi/float(N)
    for idx in range(N): 
        theta = dt*idx
        z = radius * math.cos(theta)
        y = radius * math.sin(theta)
        x = center
        print GenerateLine(x, y, z)

def SampleLine(start, end, N): 
    L = np.linalg.norm(end - start)
    dl = L/float(N-1)
    normal = (end - start) / L
    for idx in range(N): 
        point = start + normal * (dl*idx)
        print GenerateLine(point[0], point[1], point[2])

def Pool15Balls(offset, r): 
    sqrt3 = math.sqrt(3.)
    count = 1
    for row in range(4): 
        for col in range(-row, row+1, 2):
            x = offset[0] + float(row)*sqrt3*r 
            y = offset[1] + float(col)*r
            z = offset[2] + 0.0
            # print row, col, x, y, z
            ObjectPlacementTemplate(x, z, y, count)
            count += 1

if __name__ == '__main__': 
    # SampleCircle_X(0.0, 0., 100)
    SampleCircle_Z(0.0, 0.2, 10)
    # SampleLine(np.array([-0.5, 0.0, 0.0]), np.array([0.5, 0.0, 0.0]), 200)
    # Pool15Balls(np.array([0.1, 0.0, 0.0]), 0.05)
