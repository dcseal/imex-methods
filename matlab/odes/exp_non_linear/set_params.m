    % User defined parameters
    plt      = 1;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine  = 10;   % number of times to refind the grid
    mstart   = 10;   % starting resolution for the grid refinement

    tstart   = 0;
    tfinal   = 2*pi;

    meqn    = 1;

    sdc = 0;        % flag for sdc solver.  == 1 for true, o/w use rk solver
%   if( ~sdc )
        coeff_type = 'imex';
%   end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % global parameters %
    % This is an annoying way to pass parameters around %
    % I know all the functions that are being written, and nobody writes to
    % these values other than this main driving routine.

    global params

    params.sdc_order = 6;
    params.num_corrections = 2;
    params.coeffs = ['coeffs_', coeff_type];

    % problem specific parameters
    params.tfinal = tfinal;

    % might as well print the parameters used
    params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
