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

scheme_number = 2
outputdir = 'output' + names[ scheme_number ]
meqn = 2


# this variable is used for read_q()
output_num = 8

# Read in the current value of q
fname   = outputdir + ( '/q%04d.error.dat' % output_num )
fd      = open( fname, 'r' )
mt      = int( fd.readline() )

n = 0
error = np.zeros( (mt, meqn) )
for line in fd:
    line = line.split()
    for m in range(meqn):
        error[n, m] = float( line[m] )
    n = n+1
fd.flush()
fd.close()

tvec = np.linspace(0, 0.5, mt )

I = range( len(tvec) / 10 )
plt.plot( tvec[I], error[I,0], 'go' )
plt.plot( tvec[I], error[I,1], 'ro' )

plt.show()
