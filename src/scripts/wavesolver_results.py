#!/usr/bin/env python 
import numpy as np
from binary_io import *
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
        if (self.folder is None): 
            return
        filename = '%s/%s_0_listening_position.dat' %(self.folder, self.prefix)
        self.listening_points = readMatrixXdBinary(filename,1)
        self.num_listening_points = self.listening_points.shape[0]
    def Read_All_Audio(self): 
        if (self.folder is None): 
            return
        if (self.listening_points is None): 
            self.Read_Listening_Position()
        filename = '%s/%s_all_audio.dat' %(self.folder, self.prefix)
        count = 0
        with open(filename, 'rb') as stream: 
            while True:
                buf = stream.read(8*self.num_listening_points)
                if not buf: 
                    break
                try: 
                    data_step = struct.unpack('d'*self.num_listening_points, buf)
                except struct.error: 
                    # smaller than the desired size, discard this batch and break
                    break
                if count == 0: 
                    all_data = [data_step]
                else: 
                    all_data = np.append(all_data, [data_step], axis=0)
                count += 1
        return all_data
