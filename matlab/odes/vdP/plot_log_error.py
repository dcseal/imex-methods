import matplotlib

import matplotlib.pyplot as plt

import numpy as np


names = [
'_imex_eps-1.0e+00',
'_ark32_eps-1.0e+00',
'_ark43_eps-1.0e+00',
'_imex_eps-1.0e-06',
'_ark32_eps-1.0e-06',
'_ark43_eps-1.0e-06'
]

outputdir = 'output' + names[ 0 ]
numdir    = len(names)

# read in little helper file
fname = outputdir + '/qhelp.dat' 
fd = open( fname, 'r' )
nout    = int( fd.readline( ) )
meqn    = int( fd.readline()  )
tstart  = float( fd.readline() )
tfinal  = float( fd.readline() )
tau     = float( fd.readline() )
fd.flush()
fd.close()

# methods and size of epsilon
rk_type = [ 'imex', 'ark43', 'ark32' ]
eps     = [ 1e0, 1e-6 ]

er = np.zeros( (numdir, nout, meqn ) )
mt = np.zeros( nout )
dt = np.zeros( nout )

nn = 0
for na in names:

    print na
    outputdir = 'output' + na

    fname = outputdir + '/qerror.dat'
    fd    = open( fname, 'r' )

    n = 0
    for line in fd:
        line = line.split()
        mt[n] = int( line[0] )
        dt[n] = float( (tfinal - tstart) / mt[n] )
        for m in range(meqn):
            er[nn, n, m] = float( line[m+1] )
        n = n+1

    fd.close()

    nn = nn + 1

# plot the pretty pictures!
"""
plt.plot( np.log(dt), np.log( er[3,:,0] ), 'ro' )
plt.plot( np.log(dt), np.log( er[4,:,0] ), 'go' )
plt.plot( np.log(dt), np.log( er[5,:,0] ), 'bo' )
"""
plt.plot( np.log(dt), np.log( er[3,:,1] ), 'ro' )
plt.plot( np.log(dt), np.log( er[4,:,1] ), 'go' )
plt.plot( np.log(dt), np.log( er[5,:,1] ), 'bo' )

#plt.plot( np.log(dt), np.log( er[:,1] ), 'go' )

plt.show()
