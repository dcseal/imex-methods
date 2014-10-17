    % User defined parameters
    plt     = 1;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine = 10;   % number of times to refind the grid
    mxstart = 10;  % starting resolution for the grid refinement

    tstart  = 0;
    tfinal  = 0.5;

    CFL     = 0.90;     % desired cfl number


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % global parameters %
    % This is an annoying way to pass parameters around %
    % I know all the functions that are being written, and nobody writes to
    % these values other than this main driving routine.

    global params

    params.coeffs = 'coeffs_imex';

    params.meqn   = 1;  % number of equtions: this needs to be consistent with 
                        % either qinit and fluxfunc (for non-linear problems) or
                        % consistent with ConstructLinearL.

    params.maux   = 0;

    params.xlow   = 0;  % low endpoint of domain
    params.xhigh  = 1.0;  % high endpoint of domain
    %params.xhigh  = 10.0;  % high endpoint of domain

    params.sorder = 4;  % spatial order of accuracy (note: number of ghost cells
                        % comes from this)
    params.sdc_order = 4;

    % type of boundary condition
    %params.mbc_type_left  = 'periodic';
    %params.mbc_type_right = 'periodic';
    params.mbc_type_left  = 'dirichlet';
    params.mbc_type_right = 'dirichlet';
    params.mbc = 1;
    if( params.sorder > 2 )
        params.mbc = 2;
    end


    % problem specific parameters
    % params.u     = 1.0;     % advection speed
    params.u     = 0.0;     % advection speed
    params.eps   = 1e0;
    params.tfinal = tfinal;

    % might as well print the parameters used
    params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
