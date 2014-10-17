    % User defined parameters
    plt     = 1;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine = 5;   % number of times to refind the grid
    mstart  = 1;   % starting resolution for the grid refinement

    tstart  = 0;
    tfinal  = 10.0;

    meqn    = 2;

    sdc     = 0;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % global parameters %
    % This is an annoying way to pass parameters around %
    % I know all the functions that are being written, and nobody writes to
    % these values other than this main driving routine.

    global params

    params.sdc_order = 8;
    params.num_corrections = 1;
    params.coeffs = 'imex';
    params.coeffs = 'RK4';
    %params.coeffs = 'coeffs_ark2ars'


    % problem specific parameters
    params.tfinal = tfinal;

    % might as well print the parameters used
    params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
