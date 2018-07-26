#!/usr/bin/env python
import numpy as np
from binary_io import *
import glob,os
import matplotlib.pyplot as plt
################################################################################
## Class Wavesolver_Results
################################################################################
USE_CUBIC_POLY = False # use cubics to attach chunks such that they have continuous derivatives

class Wavesolver_Results:
    def __init__(self):
        self.folder = None
        self.prefix = 'test'
        # far-field listening points
        self.listening_points = None
        self.num_listening_points = 0
        # in-box sample points
        self.speaker_points = None
        self.num_speaker_points = 0
        self.box_center = None
    def Set_Folder(self, folder):
        self.folder = folder
    def Set_Prefix(self, prefix):
        self.prefix = prefix
    def Read_Listening_Position(self, TimeParallel):
        print 'Reading listening position'
        if (self.folder is None):
            return
        pattern = '%s/%s_*_listening_position.dat' %(self.folder, self.prefix)
        if (TimeParallel):
            pattern = '%s/%s_*_00000_listening_position.dat' %(self.folder, self.prefix)
        files = glob.glob(pattern)
        if not files:
            raise Exception
        filename = files[0]
        self.listening_points = readMatrixXdBinary(filename,1)
        self.num_listening_points = self.listening_points.shape[0]
        print ' %d points is read\n' %(self.num_listening_points)
    def Read_All_Audio(self, TimeParallel, NChunks, NStepsEachChunk, exclude_chunks=None, filename=None, N_speakers_along_ray=1):
        if (self.folder is None):
            return
        if (self.listening_points is None):
            self.Read_Listening_Position(TimeParallel)
        if (TimeParallel):
            all_data = np.array([])
            if exclude_chunks is not None:
                exclude_chunks = set(exclude_chunks)
            for i in range(NChunks):
                if exclude_chunks is not None and i in exclude_chunks:
                    continue
                filename = '%s/%s_%05d_all_audio.dat' %(self.folder, self.prefix, i)
                if os.path.isfile(filename):
                    filesizebytes = os.path.getsize(filename)
                else:
                    print '**WARNING** File %s does not exist' %(filename)
                    return None
                print 'Reading %s' %(filename)
                with open(filename, 'rb') as stream:
                    num_steps = int(np.floor(filesizebytes/8/self.num_listening_points))
                    if(i == 0):
                        total_steps = int(num_steps + (NChunks - 1) * NStepsEachChunk  + 19000)
                        all_data = np.zeros((total_steps, self.num_listening_points))
                    print 'Reading Chunk %d' %(i)
                    print '  reading %d steps' %(num_steps)
                    N = self.num_listening_points*num_steps
                    print 'step 0'
                    buf = stream.read(8*N)
                    print 'step 1'
                    data = np.array(struct.unpack('d'*N, buf))
                    print 'step 2'
                    data = data.reshape((num_steps, self.num_listening_points))
                    print data
                    # make chunks have continous derivatives
                    # fit a cubic spline to first and last points

                    if (USE_CUBIC_POLY):
                        d0 = data[0,:]
                        dend = data[-1,:]
                        # Set the derivatives at the endpoints to 0
                        dder0 = np.zeros((1,self.num_listening_points))
                        dderend = np.zeros((1,self.num_listening_points))
                        # fit a t^3 + bt^2 + ct + d
                        t = np.linspace(0.0,1.0, num_steps).reshape((num_steps,1))
                        t2 = np.multiply(t,t)
                        t3 = np.multiply(t2,t)
                        c = dder0
                        d = d0
                        a = c + 2 * d + dderend - 2 * dend
                        b = 3*dend - dderend - 2 * c - 3 * d
                        # multiply in the right order (tall matrix * long matrix) to fit the shape of data.
                        spline = t3 * a + t2 * b + t * c + d
                        data -= spline

                    all_data[i * NStepsEachChunk : i * NStepsEachChunk + num_steps, 0:self.num_listening_points] += data;
                    # for row in range(num_steps):
                    #     buf = stream.read(8*self.num_listening_points)
                    #     all_data[row,:] = struct.unpack('d'*self.num_listening_points, buf)
            return all_data
        else:
            if filename is None:
                filename = '%s/%s_all_audio.dat' %(self.folder, self.prefix)
            if os.path.isfile(filename):
                filesizebytes = os.path.getsize(filename)
                print 'file is %u bytes' %(filesizebytes)
            else:
                print '**WARNING** File %s does not exist' %(filename)
                return None
            print 'Reading all_audio.dat'
            with open(filename, 'rb') as stream:
                self.num_listening_points *= N_speakers_along_ray # FIXME only for delay_line_2nd
                num_steps = int(np.floor(filesizebytes/8/self.num_listening_points))
                print '  reading %d steps' %(num_steps)
                N = self.num_listening_points*num_steps
                print 'step 0'
                buf = struct.unpack('d'*N, stream.read(8*N))
                print 'step 1'
                all_data = np.array(buf)
                print 'step 2'
                all_data = all_data.reshape((num_steps, self.num_listening_points))
                print all_data
                # for row in range(num_steps):
                #     buf = stream.read(8*self.num_listening_points)
                #     all_data[row,:] = struct.unpack('d'*self.num_listening_points, buf)
            return all_data

    def Read_All_Audio_Adaptive(self, sampfreq, N_speakers_along_ray=6):
        if (self.folder is None):
            return
        if (self.listening_points is None):
            self.Read_Listening_Position(True)
        dones = sorted(glob.glob('%s/%s_*_sim-done' %(self.folder, self.prefix)))

        idxset = []
        maxtime = 0.
        for d in dones:
            idx_s = d.split('/')[-1].split('_')[1]
            f_time  = '%s/%s_%s_time_range' %(self.folder, self.prefix, idx_s)
            with open(f_time, 'r') as stream:
                lines = stream.readlines()
                maxtime = max(float(lines[1]), maxtime)
            idxset.append(idx_s)

        dt = 1./float(sampfreq)
        num_total_steps = int(np.ceil(maxtime/dt))
        self.num_listening_points *= N_speakers_along_ray
        print 'total steps = ', num_total_steps
        print 'self.num_listening_points = ', self.num_listening_points
        # tall and skinny
        all_data = np.zeros((num_total_steps, self.num_listening_points))

        for idx_s in idxset:
            print '\n\n'
            f_audio = '%s/%s_%s_all_audio.dat' %(self.folder, self.prefix, idx_s)
            f_time  = '%s/%s_%s_time_range'    %(self.folder, self.prefix, idx_s)
            assert(os.path.isfile(f_audio))
            # read time range
            with open(f_time, 'r') as stream:
                lines = stream.readlines()
                t0 = float(lines[0])
                t1 = float(lines[1])
                t2 = float(lines[2])
            # read audio data
            print 'Reading audio file: %s' %(f_audio)
            print 'Time range is ', t0, t1, t2
            filesizebytes = os.path.getsize(f_audio)
            with open(f_audio, 'rb') as stream:
                num_steps = int(np.floor(filesizebytes/8/self.num_listening_points))
                print '  reading %d steps' %(num_steps)
                N = self.num_listening_points*num_steps
                buf = stream.read(8*N)
                data = np.array(struct.unpack('d'*N, buf))
                data = data.reshape((num_steps, self.num_listening_points))
                idx = int(idx_s)

                post_process_change_termination = False
                if not post_process_change_termination:
                    row_0 = int(t0/dt)
                    if (t2 < maxtime):
                        row_1 = row_0 + data.shape[0]
                        print 'push data rows: ', row_0, row_1
                        all_data[row_0:row_1,:] = data
                    else: # last one
                        n = all_data.shape[0] - row_0
                        row_1 = row_0 + n
                        all_data[row_0:,:] = data[:n,:]
                else: # hack to post process early termination
                    newAbsLowerBound = 2E-5
                    newWindowLength = 0.025
                    if t1 == maxtime:
                        row_0 = int(t0/dt)
                        if (t2 < maxtime):
                            row_1 = row_0 + data.shape[0]
                            all_data[row_0:row_1,:] = data
                        else: # last one
                            n = all_data.shape[0] - row_0
                            row_1 = row_0 + n
                            all_data[row_0:,:] = data[:n,:]
                    else:
                        t = t0
                        count = 0
                        while t < t1:
                            count += 1
                            t += dt
                        last = t1
                        while count < num_steps:
                            if max(np.abs(data[count,:])) > newAbsLowerBound:
                                last = t
                            else:
                                if (t - last > newWindowLength):
                                    break
                            t += dt
                            count += 1
                        row_0 = int(t0/dt)
                        all_data[row_0:row_0+count,:] = data[:count,:]
        return all_data

    @staticmethod
    def ComputeFFT(dt, y):
        N = len(y)
        x = np.linspace(0.0, float(N)*dt, N)
        yf = np.fft.fft(y)[:N//2]
        xf = np.linspace(0.0, 1.0/(2.0*dt), N/2)
        magyf = np.abs(yf)
        return xf, magyf
