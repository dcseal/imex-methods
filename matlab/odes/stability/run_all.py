#!/usr/bin/env python

import os
from subprocess import Popen
from string import Template

set_params_template = Template( """
    % User defined parameters
    plt     = 0;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine = 12;   % number of times to refind the grid
    mstart  = 4;   % starting resolution for the grid refinement

    tstart  = 0;
    tfinal  = 0.50;

    coeff_type = '$coeff_type';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % global parameters %
    % This is an annoying way to pass parameters around %
    % I know all the functions that are being written, and nobody writes to
    % these values other than this main driving routine.

    global params

    params.coeffs = ['coeffs_', coeff_type];

    % problem specific parameters
    params.tfinal = tfinal;

    params.meqn   = 2;
    params.eps = $epsilon;

    % might as well print the parameters used
    params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % outputdir = '$$PWD/output_';
    outputdir = 'output_';
    outputdir = [outputdir, coeff_type, '_eps-'];
    outputdir = [outputdir, num2str( params.eps, '%1.1e' ), '/' ];
"""
)


rk_type = [
'imex',
'ark43',
'ark32'
]

eps = [
1e0,
1e-6
]

for r in rk_type:
    for e in eps:

        # open and write the set_params file
        s = set_params_template.substitute( coeff_type=r, epsilon=e )
        fd = open( 'set_params.m', 'w' )
        fd.write( s )
        fd.flush()
        fd.close()

        # run matlab with these parameters
        mycmd = 'time matlab -nodisplay -nosplash -r main'

        print 'running ' + mycmd
        p = Popen(mycmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]

"""
names = [
'_imex_eps-1.0e+00',
'_ark32_eps-1.0e+00',
'_ark43_eps-1.0e+00',
'_imex_eps-1.0e-06',
'_ark32_eps-1.0e-06',
'_ark43_eps-1.0e-06'
]

for n in names:

    print 'method:', n[1::]
    cmd = 'python run_convergence.py -o ' + 'output' + n
    os.system( cmd )
"""
