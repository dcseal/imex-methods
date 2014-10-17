#!/usr/bin/env python

import numpy as np
import scipy as sp
import math

outputdir = './output'

# read all the items in qhelp #
fd      = open( outputdir + '/qhelp.dat', 'r' )
nout    = int( fd.readline() )
meqn    = int( fd.readline() )
mx      = int( fd.readline() )
xlow    = float( fd.readline() )
xhigh   = float( fd.readline() )
dx      = float( fd.readline() )

def readfile( framenum ):

    fname = outputdir + '/q%04d.dat' % framenum
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

print 'error = %2.8e' % (sum( abs(q-q0) ) / sum( abs(q0) ) )
