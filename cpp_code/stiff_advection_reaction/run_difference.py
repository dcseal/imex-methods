import numpy as np
import scipy as sp
import math

mx = 400

def readfile( framenum ):

    fname = './output/q%04d.dat' % framenum
    fd = open( fname, 'r' )
    t = float( fd.readline() )


    q   = np.zeros( mx )
    try:
        n = 0
        for line in fd:
            q[n] = float( line )
            n    = n + 1

    finally:
        fd.close()

    return q

q0 = readfile(0)
q  = readfile(1)

print 'error = ', sum( abs(q-q0) ) / sum( abs(q) )
