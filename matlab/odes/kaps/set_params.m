
    % User defined parameters
    plt     = 1;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine = 3;   % number of times to refind the grid
    mstart  = 10;   % starting resolution for the grid refinement


    tstart  = 0;
    tfinal  = 1.0;

    sdc = 1;        % flag for sdc solver.  == 1 for true, o/w use rk solver
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
    params.num_corrections = 1;
    params.coeffs = ['coeffs_', coeff_type];

    % problem specific parameters
    params.tfinal = tfinal;

    params.meqn   = 2;
    params.eps = 1e-02;

    % might as well print the parameters used
    params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    outputdir = 'output_';
    outputdir = [outputdir, coeff_type, '_eps-'];
    outputdir = [outputdir, num2str( params.eps, '%1.1e' ), '/' ];

