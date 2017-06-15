
    % User defined parameters
    plt     = 1;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine = 10;  % number of times to refind the grid
    mstart  =  4;  % starting resolution for the grid refinement

    tstart  = 0;
    tfinal  = 4.0;

    advance = @semi_implicit_sdc_theta_integrator;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % global parameters %
    % This is an annoying way to pass parameters around %
    % I know all the functions that are being written, and nobody writes to
    % these values other than this main driving routine.

    global params

    params.sdc_order = 4;
    params.num_corrections = 3;

    % problem specific parameters
    params.tfinal = tfinal;

    params.meqn   = 2;
    params.eps    = 1.0;
%   params.eps    = 1e-03;

    params.fname = [['convergence_', func2str(advance), '_eps', num2str(abs(params.eps)), '.dat']];

    % might as well print the parameters used
    params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % outputdir = '$PWD/output_';
%   outputdir = 'output_';
%   outputdir = [outputdir, coeff_type, '_eps-'];
%   outputdir = [outputdir, num2str( params.eps, '%1.1e' ), '/' ];
