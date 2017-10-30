#! env python2

import meshio
import collections
import struct
import os
import math


def loadVtk(inFile):
    # Meshio way
    VtkData = collections.namedtuple('VtkData', 'points, cells, point_data, cell_data, field_data')
    data = VtkData._make(meshio.read(inFile))

    return data


def loadInfoFile(inFile):
    data = open(inFile).readlines()

    # skip first line (time)
    data = data[1:]

    freqs = []

    for d in data:
        parts = d.split()
        freqs.append(float(parts[1]))

    return freqs


def loadGmsh(inFile):
    inF = open(inFile).readlines()

    # skip first 4 lines
    inF = inF[4:]

    # next line is # of vertices
    numVerts = int(inF[0])
    inF = inF[1:]

    # read the vertices
    verts = []
    for i in range(numVerts):
        v = [float(x) for x in inF[i].split()[1:]]
        verts.append(v)

    inF = inF[numVerts:]

    # skip two lines
    inF = inF[2:]

    numTris = int(inF[0])
    inF = inF[1:]

    triangles = []
    triTypes = []

    for i in range(numTris):
        parts = inF[i].split()

        triTypes.append(int(parts[3]))

        triangles.append([int(x)-1 for x in parts[5:8]])

    Mesh = collections.namedtuple('Mesh', 'vertices, triangles, triTypes')
    mesh = Mesh._make((verts, triangles, triTypes))

    return mesh


def computeReordering(vtk, gmsh):
    # indexing from gmsh to vtk
    ordering = []
    for t in gmsh.triangles:
        centroid = [1./3. * sum([gmsh.vertices[x][0] for x in t]),
                    1./3. * sum([gmsh.vertices[x][1] for x in t]),
                    1./3. * sum([gmsh.vertices[x][2] for x in t]) ]

        closestDist = 1e8
        closestInd = 0

        for i in range(len(vtk.cells['triangle'])):
            vt = vtk.cells['triangle'][i]
            centroidVtk = [1./3. * sum([vtk.points[x][0] for x in vt]),
                           1./3. * sum([vtk.points[x][1] for x in vt]),
                           1./3. * sum([vtk.points[x][2] for x in vt])]

            dist = sum([(x - y)**2 for x,y in zip(centroid, centroidVtk)])

            if dist < closestDist:
                closestDist = dist
                closestInd = i

        ordering.append(closestInd)

    return ordering


if __name__ == '__main__':
    files = sorted([x for x in os.listdir('.') if x.startswith('interior') and x.endswith('vtu')])

    times = list(sorted(set([x[-14:-4] for x in files])))

    meshFiles = sorted([x for x in os.listdir('..') if x.startswith('gmsh')])
    infoFiles = sorted([x for x in os.listdir('.') if x.startswith('info')])

    RHO = 1.184

    FLUID_AIR = 1

    for t in times:
        nowFiles = sorted([x for x in files if t in x])

        meshFile = [x for x in meshFiles if t in x]

        if len(meshFile) > 1:
            raise Exception('bad length meshFile')

        meshData = loadGmsh('../' + meshFile[0])
        data = [loadVtk(f) for f in nowFiles]

        #ordering = computeReordering(data[0], meshData)
        #print(ordering)

        infoFile = [x for x in infoFiles if t in x]

        if len(meshFile) > 1:
            raise Exception('bad length info file')

        freqs = loadInfoFile(infoFile[0])

        fName = 'fullHelmSolution-' + t + '.dat'

        of = open(fName, 'wb')

        for i in range(len(data)):
            d = data[i]
            omega = 2 * math.pi * freqs[i]

            # Dummy dirichlet data
            for v in d.points:
                of.write(struct.pack('d', 0))
                of.write(struct.pack('d', 0))

            # Neumann data
            nData = d.cell_data['triangle']['n']
            for i in range(len(meshData.triangles)):
                if meshData.triTypes[i] == FLUID_AIR:
                    of.write(struct.pack('d', 0))
                    #of.write(struct.pack('d', -RHO * omega * nData[ordering[i]]))
                    of.write(struct.pack('d', -RHO * omega * nData[i]))

        of.close()
        print('finished file: ' + fName)


