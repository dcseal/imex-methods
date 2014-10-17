#!/usr/bin/env python

import os
from subprocess import Popen
from string import Template

set_params_template = Template( """
    % User defined parameters
    plt     = 0;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine = 5;   % number of times to refind the grid
    mstart  = 2;   % starting resolution for the grid refinement

    tstart  = 0;
    tfinal  = 2.0;

    sdc = $sdc_flag;        % flag for sdc solver.  == 1 for true, o/w use rk solver

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % global parameters %
    % This is an annoying way to pass parameters around %
    % I know all the functions that are being written, and nobody writes to
    % these values other than this main driving routine.

    global params

    params.sdc_order = $sdc_order;
    params.coeffs = 'coeffs_$coeff_type';

    % problem specific parameters
    params.tfinal = tfinal;

    params.meqn   = 1;
    params.tau    = $tau;

    % might as well print the parameters used
    params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    outputdir = 'output_';
    outputdir = [outputdir, '$coeff_type', '_tau-'];
    outputdir = [outputdir, num2str( params.tau, '%1.1e' ), '/' ];
"""
)


rk_type = [ 'imex', 'ark43', 'ark32' ]

tau = [ 1.0e0, 1./25.]

#for r in rk_type:
#for o in range(2,5):
for o in range(6,7):

    #r = 'sdc' + str(o)
    r = 'imex'
    for t in tau:

        # open and write the set_params file
        s = set_params_template.substitute( sdc_flag=1, sdc_order=o, coeff_type=r, tau=t )
        fd = open( 'set_params.m', 'w' )
        fd.write( s )
        fd.flush()
        fd.close()

        # run matlab with these parameters
        # mycmd = 'time matlab -nodisplay -nosplash -r main'
        mycmd = 'time octave --eval main'

        print 'running ' + mycmd
        p = Popen(mycmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
