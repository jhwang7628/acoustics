#!/usr/bin/env python
import numpy as np
import sys
################################################################################
################################################################################
class TriangleMesh:
    def __init__(self):
        self.v = []
        self.vn = []
        self.f = []
        self.misc_lines = dict()
        self.vt_lines = []
        self.f_lines = []
    def Read(self, filename):
        with open(filename, 'r') as stream:
            lines = stream.readlines()
            ii = 0
            for line in lines:
                tokens = line.split()
                if tokens[0] == 'v':
                    self.v.append([float(tokens[1]), float(tokens[2]), float(tokens[3])])
                elif tokens[0] == 'vt':
                    self.vt_lines.append(line)
                elif tokens[0] == 'vn':
                    self.vn.append([float(tokens[1]), float(tokens[2]), float(tokens[3])])
                elif tokens[0] == 'f':
                    self.f_lines.append(line)
                    self.f.append([int(tokens[1].split('/')[0])-1, \
                                   int(tokens[2].split('/')[0])-1, \
                                   int(tokens[3].split('/')[0])-1])
                else:
                    self.misc_lines[ii] = line
                ii += 1
            self.v = np.array(self.v)
            self.vn = np.array(self.vn)
            self.f = np.array(self.f, dtype=int)
    def Write(self, filename):
        with open(filename, 'w') as stream:
            count = 0
            v_written = False
            vn_written = False
            vt_written = False
            while True:
                if count in self.misc_lines:
                    stream.write(self.misc_lines[count])
                    count += 1
                    continue
                if not v_written:
                    for ii in range(self.v.shape[0]):
                        stream.write('v ')
                        for jj in range(self.v.shape[1]):
                            stream.write('%.12f ' %(self.v[ii,jj]))
                        stream.write('\n')
                        count += 1
                    v_written = True
                    continue

                if not vt_written:
                    for l in self.vt_lines:
                        stream.write(l)
                        count += 1
                    vt_written = True
                    continue

                if not vn_written:
                    for ii in range(self.vn.shape[0]):
                        stream.write('vn ')
                        for jj in range(self.vn.shape[1]):
                            stream.write('%.12f ' %(self.vn[ii,jj]))
                        stream.write('\n')
                        count +=1 
                    vn_written = True
                    continue

                # write faces and break
                f_count = 0
                while f_count < len(self.f_lines):
                    if count in self.misc_lines:
                        stream.write(self.misc_lines[count])
                    else:
                        l = self.f_lines[f_count]
                        stream.write(l)
                        f_count += 1
                    count += 1
                break

################################################################################
################################################################################
if __name__ == '__main__':
    mesh = TriangleMesh()
    mesh.Read(sys.argv[1])
    mesh.Write('test_tri.obj')
