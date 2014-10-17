    % User defined parameters
    plt     = 1;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine = 10;   % number of times to refind the grid
    mstart  = 2;   % starting resolution for the grid refinement

    tstart  = 0;
    tfinal  = 2.0;

    sdc = 1;        % flag for sdc solver.  == 1 for true, o/w use rk solver
%   if( ~sdc )
        coeff_type = 'imex';
        %coeff_type = 'ark32';
        %coeff_type = 'ark43';
%   end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % global parameters %
    % This is an annoying way to pass parameters around %
    % I know all the functions that are being written, and nobody writes to
    % these values other than this main driving routine.

    global params

    params.num_corrections = 1;
    params.sdc_order = 6;

    params.coeffs = ['coeffs_', coeff_type];

    params.tau = 0.1;

    % problem specific parameters
    params.tfinal = tfinal;

    % might as well print the parameters used
    params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
