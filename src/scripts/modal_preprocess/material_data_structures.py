#!/usr/bin/env python
class Material:
    def __init__(self):
        self.name = ''
        self.traits = dict()
    @staticmethod
    def Read_Material(mat_file):
        with open(mat_file, 'r') as stream:
            lines = stream.readlines()[-2:]
            assert(len(lines) >= 2)
            assert(lines[0][0] == '#')
            traits = lines[0].split()[1:]
            values = lines[1].split()
            assert(len(traits) == len(values))
            mat = Material()
            mat.name = '.'.join(mat_file.split('/')[-1].split('.')[:-1])
            for ii in range(len(traits)):
                mat.traits[traits[ii]] = float(values[ii])
            return mat
        return None

if __name__ == '__main__':
    glass = Material.Read_Material('/home/jui-hsien/data/models/materials/glass.txt')
    print glass.name, glass.traits
    ceramics = Material.Read_Material('/home/jui-hsien/data/models/materials/ceramics.txt')
    print ceramics.name, ceramics.traits
