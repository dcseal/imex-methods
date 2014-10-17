
    % User defined parameters
    plt     = 1;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine = 10;   % number of times to refind the grid
    mstart  = 2;   % starting resolution for the grid refinement

    tstart  = 0;
    tfinal  = 2.0;

    sdc = 0;        % flag for sdc solver.  == 1 for true, o/w use rk solver

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % global parameters %
    % This is an annoying way to pass parameters around %
    % I know all the functions that are being written, and nobody writes to
    % these values other than this main driving routine.

    global params

    params.sdc_order = 6;
    params.coeffs = 'coeffs_ark32';
    params.num_corrections = 3;

    % problem specific parameters
    params.tfinal = tfinal;

    params.meqn   = 1;
    params.tau    = 1.0e-4;

    % might as well print the parameters used
    params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    outputdir = 'output_';
    outputdir = [outputdir, 'imex', '_tau-'];
    outputdir = [outputdir, num2str( params.tau, '%1.1e' ), '/' ];
