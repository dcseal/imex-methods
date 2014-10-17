#!/usr/bin/env python

import os
from subprocess import Popen

names = [
'_imex_eps-1.0e+00',
'_ark32_eps-1.0e+00',
'_ark43_eps-1.0e+00',
'_imex_eps-1.0e-06',
'_ark32_eps-1.0e-06',
'_ark43_eps-1.0e-06'
]

for n in names:

    # latex friendly print statement here
    str = n[1::].split('_')
    meth = str[0]
    eps  = str[1].rstrip('-')
    p = str[1][ ( len( str[1] ) - 3) : len(str[1]) ]
    eps = '  $ \\eps = 10^{%d}$' % int( p )
    #print meth + eps

    #print 'method: $', n[1::], '$'
    cmd = 'python run_convergence.py -o ' + 'output' + n

    print '\\begin{figure}'
    p   = Popen( cmd, shell=True )
    sts = os.waitpid(p.pid, 0)[1]

    print '\\caption{Van der Pol test problem. ' + meth + eps + '}'
    print '\\end{figure}'
