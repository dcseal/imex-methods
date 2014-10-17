import warnings
warnings.filterwarnings("ignore")

import numpy as np
import scipy as sp
from time import time
from math import log


import os
import sys
from optparse import OptionParser

def main(outputdir):

    fname = outputdir + '/qhelp.dat' 
    fd = open( fname, 'r' )
    nout    = int( fd.readline( ) )
    meqn    = int( fd.readline()  )
    tstart  = float( fd.readline() )
    tfinal  = float( fd.readline() )
    tau     = float( fd.readline() )
    fd.flush()
    fd.close()


    def write_pointwise_err( n, q, qex ):

        fname = outputdir + ( '/q%04d.error.dat' % (n+1) )
        fd = open( fname, 'w' )
        fd.write( str( len(q) ) + '\n' )
        for n in range( len( q ) ):
            for m in range(meqn):
                val = '%2.15e ' % ( q[n,m] - qex[n,m] )
                fd.write( val )
            fd.write('\n')
        fd.flush()
        fd.close()



    # this variable is used for read_q()
    output_num = 1
    def read_q():

        # Read in the current value of q
        fname = outputdir + ( '/q%04d.dat' % output_num )
        fd = open( fname, 'r' )
        mt = int( fd.readline() )

        n = 0
        q = np.zeros( (mt, meqn) )
        for line in fd:
            line = line.split()
            for m in range(meqn):
                q[n, m] = float( line[m] )
            n = n+1
        fd.flush()
        fd.close()

        return q, mt

    # read in a bunch of earlier values of q
    q      = []
    mt_vec = []
    output_num = 1
    for n in range(nout):

        qn, mtn = read_q()
        q.append( qn )
        mt_vec.append( mtn )
        output_num = output_num + 1
        mtmax = mtn


    # read in the most refined value of q ('qex'):
    dt   = (tfinal-tstart) / mtmax
    tvec = dt * np.arange(mtmax)
    qex_big      = np.zeros( (mtmax,2) )
    qex_big[:,0] = 1/(1+tau**2)*( -np.exp(-tvec/tau) + tau*np.sin(tvec) + np.cos(tvec) )
    
    # print all the errors 
    numerr = int(nout)
    er = np.zeros( (numerr, 2) )


    fname = outputdir + '/qerror.dat'
    fd = open( fname, 'w' )

    print ' '
    print '\\[ \\begin{array}{c|c|c}'
    print '    mt & error( q ) &  log2( ratio ) \\\\'
    print '\\hline'

    for n in range( numerr ):

        qn   = q[n]
        mt   = mt_vec[n]
        skip = mtmax / mt

        qex = qex_big[::skip]
        
        write_pointwise_err( n, qn, qex )

        # vdP problem shows order reduction in one variable only
        er[n,0] = np.linalg.norm( qex[:,0] - qn[:,0], 1 ) / np.linalg.norm( qex[:,0], 1 )
        if( n > 0 ):
            l0 = log( er[n-1,0] / er[n,0], 2)
#           print '    log2(ratio0) = %2.5f' % (l0) + \
#           '    mt = %4d; er[0] = %2.3e' %( mt, er[n,0])

            print '    %4d & %2.3e ' % (mt, er[n,0] ) + \
            '    & %2.3f \\\\' % ( l0 )

        # write total error to file
        val = '%d %2.15e %2.15e\n' % ( mt, er[n,0], er[n,1] )
        fd.write( val )


    fd.flush()
    fd.close()

    print '\\end{array} \\]'

if __name__ == "__main__":
    """Write a docstring here.  No seriously, you should write a docstring!"""

    parser = OptionParser()
    parser.add_option("-o", "--outputdir", dest="outputdir",
                      help="output directory to be created",
                      default='output', metavar="OUTPUTDIR")

    (options, args) = parser.parse_args()


    main( options.outputdir )
