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

framenum = 0
outputdir = 'output' + names[ framenum ]

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


# this variable is used for read_q()
output_num = 8

# Read in the current value of q
fname = outputdir + ( '/q%04d.error.dat' % output_num )
fd = open( fname, 'r' )
mt = int( fd.readline() )

n = 0
error = np.zeros( (mt, meqn) )
for line in fd:
    line = line.split()
    for m in range(meqn):
        error[n, m] = float( line[m] )
    n = n+1
fd.flush()
fd.close()

dt = (tfinal-tstart) / mt
tvec = dt * np.arange(mt)

# plot the pretty pictures!
plt.plot( tvec, error[:,0], 'go' )
plt.plot( tvec, error[:,1], 'ro' )

plt.show()
