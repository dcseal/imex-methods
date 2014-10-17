%% parameters file %%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % local parameters %
    % These are only used in the main driver
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % User defined parameters
    plt     = 1;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine = 1;      % number of times to refind the grid
    mxstart = 50;     % where to start the refinement

    tstart  = 0;
    tfinal  = 0.70;

    CFL     = 0.9;     % desired cfl number


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % global parameters %
    % This is an annoying way to pass parameters around %
    % I know all the functions that are being written, and nobody writes to
    % these values other than this main driving routine.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global params

    params.coeffs = 'coeffs_imex';
    params.sdc    = 1;

    params.meqn   = 2;  % number of equtions: this needs to be consistent with 
                        % either qinit and fluxfunc (for non-linear problems) or
                        % consistent with ConstructLinearL.

    params.maux   = 1;  % number of auxilary variables

    params.xlow   = 0;  % low endpoint of domain
    params.xhigh  = 1;  % high endpoint of domain

    params.sorder = 6;  % spatial order of accuracy (note: number of ghost cells
                        % comes from this)
    params.sdc_order = 6;  % spatial order of accuracy (note: number of ghost cells

    % type of boundary condition
    params.mbc_type_left  = 'periodic';
    params.mbc_type_right = 'periodic';

    % number of boundary cells depends on sorder only.
    params.mbc = 1;
    if( params.sorder > 2 )
        params.mbc = 3;
    end


    % problem specific parameters
    params.u     = 1.0;     % advection speed
    params.tau   = 1e-2;
    params.tf    = tfinal;

    % might as well print the parameters used
    params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
