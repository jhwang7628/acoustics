#!/usr/bin/env python
import numpy as np
import struct
import scipy.io.wavfile as wavfile

################################################################################
################################################################################
def readVectorBinary(filename):
    ifile = open( filename, 'rb' )
    N1 = np.fromfile( ifile, dtype=np.int32, count=1 )
    print('read vector size: (%u)' %(N1))
    tmp = np.fromfile( ifile, dtype=float, count=N1 )
    return tmp

################################################################################
################################################################################
def readMatrixXdBinary(filename, verbose=1):
    ifile = open( filename, 'rb' )
    N1 = np.fromfile( ifile, dtype=np.int32, count=1 )[0]
    N2 = np.fromfile( ifile, dtype=np.int32, count=1 )[0]
    if (verbose >= 1):
        print('read matrix size: (%u x %u)' %(N1,N2))
    tmp = np.fromfile( ifile, dtype=float, count=N1*N2 )
    tmp = tmp.reshape( (N1, N2) )
    return tmp

################################################################################
################################################################################
def writeVectorBinary(filename, A):
    print('write vector size: (%u)' %(A.shape))
    ofile = open( filename, 'wb' )
    if ofile:
        ofile.write( struct.pack('i', int(A.shape[0])) )
        A.tofile( ofile )
        ofile.close()
    else:
        print('error')

################################################################################
################################################################################
def writeMatrixXdBinary(filename, A):
    print('write matrix size: (%u x %u)' %(A.shape[0],A.shape[1]))
    ofile = open( filename, 'wb' )
    if ofile:
        ofile.write( struct.pack('i', int(A.shape[0])) )
        ofile.write( struct.pack('i', int(A.shape[1])) )
        A.tofile( ofile )
        ofile.close()
    else:
        print('error')

################################################################################
################################################################################
if __name__ == "__main__":
    A = np.random.randn( 10,10 )
    dat='tmp.dat'
    writeMatrixXdBinary( dat, A.copy() )
    AA = readMatrixXdBinary( dat )
