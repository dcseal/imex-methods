    % User defined parameters
    plt     = 1;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine = 1;      % number of times to refind the grid

    mxstart = 100;
    mystart = 100;


    tstart  = 0;
    tfinal  = 1.0;

    CFL     = 0.95;     % desired cfl number


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % global parameters %
    % This is an annoying way to pass parameters around %
    % I know all the functions that are being written, and nobody writes to
    % these values other than this main driving routine.

    global params

    params.coeffs = 'coeffs_imex';
    params.coeffs = 'coeffs_rk';
    %params.coeffs = 'coeffs_ark43';
    %params.coeffs = 'coeffs_ark32';

    params.meqn   = 1;  % number of equtions: this needs to be consistent with 
                        % either qinit and fluxfunc (for non-linear problems) or
                        % consistent with ConstructLinearL.

    params.maux   = 0;

    params.xlow   = 0;  % low endpoint of domain
    params.xhigh  = 1;  % high endpoint of domain

    params.ylow   = 0;  % low endpoint of domain
    params.yhigh  = 1;  % high endpoint of domain

    params.sorder = 4;  % spatial order of accuracy (note: number of ghost cells
                        % comes from this)

    % number of boundary cells depends on sorder only.
    params.mbc_type_left  = 'periodic';
    params.mbc_type_right = 'periodic';

    params.mbc = 1;
    if( params.sorder > 2 )
        params.mbc = 2;
    end


    % problem specific parameters
    params.u     =  1.0;     % advection speed
    params.v     = -1.0;     % advection speed
    params.tfinal = tfinal;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
