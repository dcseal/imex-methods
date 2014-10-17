#!/usr/bin/env python
from subprocess import call
import numpy as np
import scipy as sp
from scipy.integrate import quad
import math
import sys

def qexact( x ):
    return 1.0 + np.sin(10.0*np.pi*x)

nout = 10
error = np.zeros(nout)
for n in range(nout):

    mx  = 10 * 2**n;
    node = np.linspace(0,1,mx+1)
    dx  = node[1] - node[0]

    q = np.zeros(mx)
    for i in range(mx):
        q[i] = quad( qexact, node[i], node[i+1], epsabs=1e-15 )[0] / dx

    # remove the output file
    call( 'rm out.dat', shell=True )

    # read in the output file
    cmd = './main %d > out.dat' % mx
    call( cmd, shell=True );

    fd = open( 'out.dat', 'r' )
    er = 0
    for i in range(mx):
        line = float( fd.readline() )

    #   str = 'my quadrature = %2.5e, exact answer = %2.5e'  % (line,q[i])
    #   str = str + '   adding in error = %2.5e' % abs(line-q[i])
    #   print str

        er = er + abs( line - q[i] )
        
    error[n] = er / sum( abs(q) )

    if( n > 0 ):
        sys.stdout.write( 'error[n] = %2.8e' % ( error[n] ) )
        print '    log2(ratios) = %2.5f' % math.log( error[n-1] / error[n], 2 )
