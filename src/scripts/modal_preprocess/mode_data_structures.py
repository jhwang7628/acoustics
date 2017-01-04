#!/usr/bin/env python
import os,struct,math
import numpy as np
################################################################################
################################################################################
class Modes:
    def __init__(self): 
        self.num_modes = 0
        self.num_vertices = 0
        self.stride = 0    # stride = 3*number of vertices
        self.eigenvalues = None
        self.eigenvectors = None

    def Get_All_Frequencies(self, density):
        frequencies = np.zeros(self.num_modes)
        for i in range(self.num_modes):
            frequencies[i] = np.sqrt(self.eigenvalues[i]/density)/(2.0*np.pi)
        return frequencies

    def Has_Modes(self):
        return self.num_modes>0

    def Write_Modes(self, out_mode_file):
        ofs = open(out_mode_file,'wb')
        ofs.write(struct.pack('i', self.stride))
        ofs.write(struct.pack('i', self.num_modes))
        for m in range(self.num_modes):
            ofs.write(struct.pack('d', self.eigenvalues[m]))
        count_inner = 0
        for iv in range(self.num_vertices):
            for m in range(self.num_modes):
                for d in range(3): 
                    ofs.write(struct.pack('d', self.eigenvectors[count_inner,d]))
                count_inner += 1
        ofs.close()
        print 'Wrote %d modes in %d stride to file: %s' %(self.num_modes, self.stride, out_mode_file)

    def Remove_Leading_Modes(self, num_leading_modes_to_remove):
        assert(num_leading_modes_to_remove <= self.num_modes)
        new_num_modes = self.num_modes - num_leading_modes_to_remove
        self.eigenvalues = self.eigenvalues[num_leading_modes_to_remove:]
        new_eigenvectors = np.zeros((self.num_vertices*new_num_modes, 3))
        new_eigenvectors[:,:] = self.eigenvectors[num_leading_modes_to_remove*self.num_vertices:,:]
        self.eigenvectors = new_eigenvectors
        self.num_modes = new_num_modes
        print '%d leading modes removed, new object has %d modes' %(num_leading_modes_to_remove, self.num_modes)

    ############################################################################
    @staticmethod
    def Read_Modes(mode_file, read_eigenvectors=True):
        assert(os.path.isfile(mode_file))
        modes = Modes()
        ifs = open(mode_file, 'rb')
        modes.stride = struct.unpack('i', ifs.read(4))[0]
        modes.num_modes = struct.unpack('i', ifs.read(4))[0] 
        modes.num_vertices = modes.stride / 3
        modes.eigenvalues = np.zeros(modes.num_modes)
        modes.eigenvectors = np.zeros((modes.num_vertices*modes.num_modes, 3))
        for m in range(modes.num_modes):
            modes.eigenvalues[m] = struct.unpack('d', ifs.read(8))[0]
        if read_eigenvectors:
            count_inner = 0
            for iv in range(modes.num_vertices):
                for m in range(modes.num_modes):
                    for d in range(3): 
                        modes.eigenvectors[count_inner,d] = struct.unpack('d', ifs.read(8))[0]
                    count_inner += 1
            ifs.close()
        print 'Read %d modes in %d stride from file: %s' %(modes.num_modes, modes.stride, mode_file)
        return modes

    @staticmethod
    def Print_Eigen_Values_From_File(modes_file, density): 
        ifs = open(modes_file,'rb')
        nsz = struct.unpack('i', ifs.read(4))[0]
        nev = struct.unpack('i', ifs.read(4))[0]
        print '#eigenvalues = %d' %(nev)
        for e in range(nev): 
            omega_squared_times_density = struct.unpack('d', ifs.read(8))[0]
            frequency = math.sqrt(omega_squared_times_density/density)/(2.0*math.pi)
            print 'Mode %d has frequency: %f Hz' %(e, frequency)
