#!/usr/bin/env python

import os
from subprocess import Popen

names = [
'_imex_tau-1.0e+00',
'_ark32_tau-1.0e+00',
'_ark43_tau-1.0e+00',
'_imex_tau-4.0e-02',
'_ark32_tau-4.0e-02',
'_ark43_tau-4.0e-02'
]

names = [
'_sdc2_tau-1.0e+00',
'_sdc3_tau-1.0e+00',
'_sdc4_tau-1.0e+00',
'_sdc2_tau-4.0e-02',
'_sdc3_tau-4.0e-02',
'_sdc4_tau-4.0e-02'
]

for n in names:

    # latex friendly print statement here
    str = n[1::].split('_')
    meth = str[0]
    tau  = str[1].rstrip('-')
    p = str[1][ ( len( str[1] ) - 3) : len(str[1]) ]
    tau = '  $ \\tau = 10^{%d}$' % int( p )

    #print 'method: $', n[1::], '$'
    cmd = 'python run_convergence.py -o ' + 'output' + n


    print '\\begin{figure}'
    p   = Popen( cmd, shell=True )
    sts = os.waitpid(p.pid, 0)[1]

    print '\\caption{Convergence table for ODE test problem. ' + meth + tau + '}'
    print '\\end{figure}'
