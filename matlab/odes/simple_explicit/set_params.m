    % User defined parameters
    plt     = 1;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine = 10;   % number of times to refind the grid
    mstart  = 20;   % starting resolution for the grid refinement

    tstart  = 0;
    tfinal  = 10.0;

    meqn    = 1;

    sdc     = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % global parameters %
    % This is an annoying way to pass parameters around %
    % I know all the functions that are being written, and nobody writes to
    % these values other than this main driving routine.

    global params

    params.sdc_order = 5;
    params.num_corrections = 4;
    params.coeffs = 'imex';
    %params.coeffs = 'coeffs_ark2ars'

    % problem specific parameters
    params.tfinal = tfinal;
    params.lambda = -10.0;

    % might as well print the parameters used
    params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
