#!/usr/bin/env python 
import numpy as np
from binary_io import *
import glob,os
################################################################################
## Class Wavesolver_Results
################################################################################
class Wavesolver_Results: 
    def __init__(self):
        self.folder = None
        self.prefix = 'test'
        self.listening_points = None
        self.num_listening_points = 0
    def Set_Folder(self, folder): 
        self.folder = folder
    def Set_Prefix(self, prefix): 
        self.prefix = prefix
    def Read_Listening_Position(self): 
        print 'Reading listening position'
        if (self.folder is None): 
            return
        pattern = '%s/%s_*_listening_position.dat' %(self.folder, self.prefix)
        files = glob.glob(pattern)
        if not files: 
            raise Exception
        filename = files[0]
        self.listening_points = readMatrixXdBinary(filename,1)
        self.num_listening_points = self.listening_points.shape[0]
        print ' %d points is read\n' %(self.num_listening_points)
    def Read_All_Audio(self): 
        if (self.folder is None): 
            return
        if (self.listening_points is None): 
            self.Read_Listening_Position()
        filename = '%s/%s_all_audio.dat' %(self.folder, self.prefix)
        if os.path.isfile(filename):
            filesizebytes = os.path.getsize(filename)
        else: 
            print '**WARNING** File %s does not exist' %(filename)
            return None
        print 'Reading all_audio.dat'
        with open(filename, 'rb') as stream: 
            num_steps = int(np.floor(filesizebytes/8/self.num_listening_points))
            print '  reading %d steps' %(num_steps)
            N = self.num_listening_points*num_steps
            print 'step 0'
            buf = stream.read(8*N)
            print 'step 1'
            all_data = np.array(struct.unpack('d'*N, buf))
            print 'step 2'
            all_data = all_data.reshape((num_steps, self.num_listening_points))
            print all_data
            # for row in range(num_steps): 
            #     buf = stream.read(8*self.num_listening_points)
            #     all_data[row,:] = struct.unpack('d'*self.num_listening_points, buf)
        return all_data

    @staticmethod
    def ComputeFFT(dt, y): 
        N = len(y)
        x = np.linspace(0.0, float(N)*dt, N)
        yf = np.fft.fft(y)[:N//2]
        xf = np.linspace(0.0, 1.0/(2.0*dt), N/2)
        magyf = np.abs(yf)
        return xf, magyf
