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
        pattern = '%s/%s_listening_position.dat' %(self.folder, self.prefix)
        files = glob.glob(pattern)
        if not files: 
            raise Exception
        filename = files[0]
        self.listening_points = readMatrixXdBinary(filename,1)
        self.num_listening_points = self.listening_points.shape[0]
    def Read_All_Audio(self): 
        print 'Reading all_audio.dat'
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
        with open(filename, 'rb') as stream: 
            num_steps = int(np.floor(filesizebytes/8/self.num_listening_points))
            all_data = np.zeros((num_steps, self.num_listening_points))
            for row in range(num_steps): 
                buf = stream.read(8*self.num_listening_points)
                all_data[row,:] = struct.unpack('d'*self.num_listening_points, buf)
        return all_data
