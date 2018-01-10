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

    def Read_All_Audio_Adaptive(self, sampfreq):
        if (self.folder is None):
            return
        if (self.listening_points is None):
            self.Read_Listening_Position(True)
        dones = glob.glob('%s/%s_*_sim-done' %(self.folder, self.prefix))

        idxset = set()
        maxtime = 0.
        for d in dones:
            idx_s = d.split('/')[-1].split('_')[1]
            f_time  = '%s/%s_%s_time_range' %(self.folder, self.prefix, idx_s)
            with open(f_time, 'r') as stream:
                lines = stream.readlines()
                maxtime = max(float(lines[1]), maxtime)
            idxset.add(idx_s)

        dt = 1./float(sampfreq)
        num_total_steps = int(np.ceil(maxtime/dt))
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

    def Read_Listening_Position_Delay_Line(self, chunkId_s):
        assert(self.folder is not None)
        filename = glob.glob('%s/%s_*_%s_listening_model_metadata' %(self.folder, self.prefix, chunkId_s))[0]
        lines = open(filename,'r').readlines()

        ii = 0
        while ii < len(lines):
            if lines[ii][0] != '#':
                break
            ii += 1
        l = lines[ii]
        self.num_listening_points = int(l.split()[0])
        self.listening_points = np.zeros((self.num_listening_points,3))
        ii += 1
        for jj in range(self.num_listening_points):
            l = lines[ii]
            tokens = l.split()
            for kk in range(3):
                self.listening_points[jj,kk] = float(tokens[kk])
            ii += 1
        self.num_speaker_points = int(lines[ii].split()[0])
        self.speaker_points = np.zeros((self.num_speaker_points,3))
        ii += 1
        for jj in range(self.num_speaker_points):
            l = lines[ii]
            tokens = l.split()
            for kk in range(3):
                self.speaker_points[jj,kk] = float(tokens[kk])
            ii += 1
        N_box_center = int(lines[ii].split()[0])
        ii += 1
        self.box_center = np.array([float(lines[ii].split()[0]),
                                    float(lines[ii].split()[1]),
                                    float(lines[ii].split()[2])])
        print 'box_center =', self.box_center

    def Read_All_Audio_Adaptive_Delay_Line_2nd(self, sampfreq):
        if (self.folder is None):
            return
        self.Read_Listening_Position_Delay_Line('00000')
        # FIXME debug START
        # dones = glob.glob('%s/%s_*_sim-done' %(self.folder, self.prefix))
        dones = glob.glob('%s/%s_*_all_audio.dat' %(self.folder, self.prefix))
        # FIXME debug END

        idxset = set()
        maxtime = 7.0 # FIXME debug
        for d in dones:
            idx_s = d.split('/')[-1].split('_')[1]
            f_time  = '%s/%s_%s_time_range' %(self.folder, self.prefix, idx_s)
            if os.path.isfile(f_time):
                with open(f_time, 'r') as stream:
                    lines = stream.readlines()
                    maxtime = max(float(lines[1]), maxtime)
            idxset.add(idx_s)

        dt = 1./float(sampfreq)
        num_total_steps = int(np.ceil(maxtime/dt))
        # tall and skinny
        all_data = np.zeros((num_total_steps, self.num_listening_points))

        for idx_s in idxset:
            print '\n\n'
            f_audio = '%s/%s_%s_all_audio.dat' %(self.folder, self.prefix, idx_s)
            f_time  = '%s/%s_%s_time_range'    %(self.folder, self.prefix, idx_s)
            assert(os.path.isfile(f_audio))
            # read time range
            if os.path.isfile(f_time):
                with open(f_time, 'r') as stream:
                    lines = stream.readlines()
                    t0 = float(lines[0])
                    t1 = float(lines[1])
                    t2 = float(lines[2])
            else:
                f_time = '%s/%s_%s_start_time' %(self.folder, self.prefix, idx_s)
                with open(f_time, 'r') as stream:
                    lines = stream.readlines()
                    t0 = float(lines[0])
                    t2 = maxtime - 1E-6
            # read audio data
            print 'Reading audio file: %s' %(f_audio)
            filesizebytes = os.path.getsize(f_audio)
            with open(f_audio, 'rb') as stream:
                idx = int(idx_s)
                self.Read_Listening_Position_Delay_Line(idx_s)
                num_steps = int(np.floor(filesizebytes/8/self.num_speaker_points))
                print '  reading %d steps' %(num_steps)
                N = self.num_speaker_points*num_steps
                buf = stream.read(8*N)
                data = np.array(struct.unpack('d'*N, buf))
                data = data.reshape((num_steps, self.num_speaker_points))

                # read box center
                x_c = self.box_center
                c = 343.0

                print 'Reconstructing listener pressure'
                for ii in range(self.num_listening_points):
                    x_L0 = self.listening_points[ii]
                    x_0 = self.speaker_points[ii*2  ]
                    x_1 = self.speaker_points[ii*2+1]

                    alpha = round(np.linalg.norm(x_L0 - x_c)/c/dt)
                    # jiggle x_L a bit so its interger multiple of stepsize
                    x_L = x_c + (x_L0-x_c)/np.linalg.norm(x_L0-x_c)*alpha*c*dt

                    print '  x_L =', x_L
                    print '  x_0 =', x_0
                    print '  x_1 =', x_1
                    r_L = np.linalg.norm(x_L - x_c)
                    r_0 = np.linalg.norm(x_0 - x_c)
                    r_1 = np.linalg.norm(x_1 - x_c)
                    dt0 = int((r_L - r_0)/c/dt)
                    dt1 = int((r_L - r_1)/c/dt)
                    print dt0, dt1
                    p_0 = data[:,ii*2  ]
                    p_1 = data[:,ii*2+1]

                    # plt.figure()
                    # plt.plot(p_0, label='p0')
                    # plt.plot(p_1, label='p1')
                    # plt.show()

                    p_L = np.zeros(num_steps)
                    dr = r_0 - r_1
                    s_0 =  r_0**2/(r_L*dr) - r_0**2*r_1/(r_L**2*dr)
                    s_1 = -r_1**2/(r_L*dr) + r_0*r_1**2/(r_L**2*dr)
                    print s_0, s_1, s_0+s_1
                    # p_L = s_0*p_0 + s_1*p_1
                    p_L[dt0:num_steps]  = s_0*p_0[:num_steps-dt0]
                    # p_L[dt1:num_steps] += s_1*p_1[:num_steps-dt1]
                    row_0 = int(t0/dt)
                    all_data[row_0:row_0+len(p_L),ii] += p_L

                # FIXME debug START
                # if (t2 < maxtime):
                #     row_1 = row_0 + data.shape[0]
                #     all_data[row_0:row_1,:] = data
                # else: # last one
                #     n = all_data.shape[0] - row_0
                #     row_1 = row_0 + n
                #     all_data[row_0:,:] = data[:n,:]
                # FIXME debug END
        return all_data

    @staticmethod
    def ComputeFFT(dt, y):
        N = len(y)
        x = np.linspace(0.0, float(N)*dt, N)
        yf = np.fft.fft(y)[:N//2]
        xf = np.linspace(0.0, 1.0/(2.0*dt), N/2)
        magyf = np.abs(yf)
        return xf, magyf
