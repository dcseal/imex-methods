
    % User defined parameters
    plt     = 1;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine = 10;   % number of times to refind the grid
    mstart  = 10;   % starting resolution for the grid refinement

    tstart  = 0;
    tfinal  = pi;

    % make sure this stuff is consistent with the stuff present in main.m
    sdc = 1;        % flag for sdc solver.  == 1 for true, o/w use rk solver
    coeff_type = 'imex';


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % global parameters %
    % This is an annoying way to pass parameters around %
    % I know all the functions that are being written, and nobody writes to
    % these values other than this main driving routine.

    global params

    params.sdc_order = 5;
    params.num_corrections = 4;
    params.coeffs = ['coeffs_', coeff_type];

    % problem specific parameters
    params.tfinal = tfinal;

    params.meqn   = 2;

    % might as well print the parameters used
    params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
