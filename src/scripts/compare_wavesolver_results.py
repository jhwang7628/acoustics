#!/usr/bin/env python
import sys,os
import matplotlib.pyplot as plt
from wavesolver_results import *
from ConfigParser import SafeConfigParser,NoOptionError
import scipy
from scipy import signal

if len(sys.argv) != 3: 
    print '**Usage: %s <config_ini_file> <config_ini_file>' %sys.argv[0]
    sys.exit()
out_config = sys.argv[1]
if not os.path.isfile(out_config): 
    print '**ERROR** Config file %s does not exist' %(out_config)
    sys.exit()

out2_config = sys.argv[2]
if not os.path.isfile(out2_config): 
    print '**ERROR** Config file %s does not exist' %(out2_config)
    sys.exit()
##
def ReqParse(parser, section, tag, val_type = 's'): 
    if val_type == 'i':   # int
        tmp = parser.getint(section, tag)
    elif val_type == 'f': # float
        tmp = parser.getfloat(section, tag)
    elif val_type == 'b': # boolean
        tmp = parser.getboolean(section, tag)
    elif val_type == 's': # string
        tmp = parser.get(section, tag)
    elif val_type == 'li': # list of int (comma separated)
        tmp = map(int, parser.get(section, tag).split(','))
    else: 
        raise TypeError
    print 'Parsing \'%s\' in \'%s\': ' %(tag, section), tmp
    return tmp

def OptParse(parser, section, tag, val_type = 's', val_default = None): 
    try: val = ReqParse(parser, section, tag, val_type)
    except NoOptionError: val = val_default
    return val
##
parser = SafeConfigParser()
parser.read(out_config)

data_dir = ReqParse(parser, 'general', 'data_dir', 's')
print data_dir
out_dir = ReqParse(parser, 'general', 'out_dir', 's')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

time_parallel = OptParse(parser, 'general', 'time_parallel', 'b', False)
results = Wavesolver_Results()
results.Set_Folder(data_dir)
all_data = np.array([])
NChunks = 1
NStepsEachChunk = 0
if(time_parallel):
    NChunks = ReqParse(parser, 'general', 'time_chunks', 'i')
    NStepsEachChunk = ReqParse(parser, 'general', 'time_steps_each_chunk', 'i')
all_data = results.Read_All_Audio(time_parallel, NChunks, NStepsEachChunk)


parser = SafeConfigParser()
parser.read(out2_config)

data_dir = ReqParse(parser, 'general', 'data_dir', 's')
print data_dir
out_dir = ReqParse(parser, 'general', 'out_dir', 's')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

time_parallel = OptParse(parser, 'general', 'time_parallel', 'b', False)
results = Wavesolver_Results()
results.Set_Folder(data_dir)
all_data2 = np.array([])
NChunks = 1
NStepsEachChunk = 0
if(time_parallel):
    NChunks = ReqParse(parser, 'general', 'time_chunks', 'i')
    NStepsEachChunk = ReqParse(parser, 'general', 'time_steps_each_chunk', 'i')
all_data2 = results.Read_All_Audio(time_parallel, NChunks, NStepsEachChunk)


##

N_points = all_data.shape[1]
N_steps  = all_data.shape[0] 

N_points2 = all_data2.shape[1]
N_steps2  = all_data2.shape[0] 
sampfreq = ReqParse(parser, 'general', 'sampfreq', 'i')
if (N_steps2 < N_steps):
    N_steps = N_steps2
all_data_diff = all_data[0:N_steps,0:N_points] - all_data2[0:N_steps,0:N_points2]
t = np.array(range(N_steps))#/float(sampfreq)
plt.figure()
plt.plot(t, all_data_diff[:,0])


plt.figure()
plt.plot(t, all_data[0:N_steps,0])
plt.plot(t, all_data2[0:N_steps,0])


plt.figure()
diffsum = np.sum( np.abs( all_data_diff), axis=None)
datasum = np.sum( np.abs( all_data[0:N_steps, :]), axis=None)
print diffsum / datasum
#plt.plot(t, np.divide( all_data_diff[:, 0] , all_data[0:N_steps,0] ) )
print " max difference: " + str(np.max(np.abs(all_data_diff), axis=None))


if ReqParse(parser, 'general', 'plot', 'b'):
    print '\n------ PLOTTING ------'
    xaxis_frame = ReqParse(parser, 'plot', 'xaxis_frame', 'b')
    xmax        = OptParse(parser, 'plot', 'xmax'       , 'f')
    plot_points = OptParse(parser, 'plot', 'plot_points', 'li')
    plt.figure()
    if (plot_points is None): ii_range = range(N_points)
    else                    : ii_range = plot_points
    for ii in ii_range: 
        if xaxis_frame: t = np.array(range(len(all_data[:,ii])))
        else:           t = np.array(range(len(all_data[:,ii])))/float(sampfreq)
        y = all_data[:,ii]
        if xmax is not None: 
            cut = -1
            for ii in range(len(t)): 
                if t[ii] > xmax: 
                    cut = ii
                    break
            if cut >= 0: 
                t = t[:cut]
                y = y[:cut]
        plt.plot(t, y)
    plt.show()

if ReqParse(parser, 'general', 'write_wav', 'b'):
    print '\n------ WRITING WAV FILES ------'
    wavfreq  = ReqParse(parser, 'wav'    , 'wavfreq' , 'i')
    prefix   = ReqParse(parser, 'wav'    , 'prefix'  , 's')
    rateRatio = float(sampfreq) / float(wavfreq)
    for ii in range(N_points): 
        print 'point %u' %(ii)
        outputdata = signal.resample(all_data[:,ii], int(float(N_steps)/rateRatio))
        normalization = np.absolute(outputdata).max()
        scipy.io.wavfile.write('%s/%s_%u.wav' %(out_dir, prefix, ii), wavfreq, outputdata/normalization)
