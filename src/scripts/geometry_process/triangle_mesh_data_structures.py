#!/usr/bin/env python
import numpy as np
class Triangle_Mesh:
    def __init__(self):
        self.V = []
        self.F = []
    @staticmethod
    def Read_Obj(objfile):
        mesh = Triangle_Mesh()
        with open(objfile, 'r') as stream:
            lines = stream.readlines()
            for line in lines:
                tokens = line.split()
                if tokens[0] == 'v':
                    mesh.V.append([float(x) for x in tokens[1:4]])
                elif tokens[0] == 'f':
                    mesh.F.append([int(x)-1 for x in tokens[1:4]])
        return mesh

if __name__ == '__main__':
    mesh = Triangle_Mesh.Read_Obj('/home/jui-hsien/data/models/A/A.obj')
    print np.array(mesh.V)
    print np.array(mesh.F, dtype=int)
