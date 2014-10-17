    clear all

    % User defined parameters
    plt     = 1;   % whether or not to turn on the plotter.  set to non-zero if
                   % yes

    nrefine = 10;      % number of times to refind the grid

    tstart  = 0;
    tfinal  = 1.0;

    CFL     = 9.0;     % desired cfl number


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % global parameters %
    % This is an annoying way to pass parameters around %
    % I know all the functions that are being written, and nobody writes to
    % these values other than this main driving routine.

    global params

    params.meqn   = 1;  % number of equtions: this needs to be consistent with 
                        % either qinit and fluxfunc (for non-linear problems) or
                        % consistent with ConstructLinearL.

    params.maux   = 0;

    params.xlow  = 0;   % low endpoint of domain
    params.xhigh = 1;   % high endpoint of domain

    params.sorder = 4;  % spatial order of accuracy (note: number of ghost cells
                        % comes from this)

    % number of boundary cells depends on sorder only.
    params.mbc_type = 'periodic';
    params.mbc = 1;
    if( params.sorder > 2 )
        params.mbc = 2;
    end


    % problem specific parameters
    params.u     = 1.0;     % advection speed

    % split problem:  here (uimp + uexp) = u
    params.uimp  = 0.9;
    params.uexp  = 0.1;

%   params.r     = 0.2;     % parameter in one of the test problems
%   params.tau   = 1e-9;

    % might as well print the parameters used
    params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Outer loop.  Each time this halves the grid size if nrefine > 1 %
    er = zeros(nrefine,1);
    for no=1:nrefine


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create new grid
        %
        mx   = 10*2^(no-1);  % number of grid cells (not counting boundary cells)
        params.mx = mx;


        % set up grid spacing and nodal points
        xhigh = params.xhigh;  xlow = params.xlow;  mbc = params.mbc;
        dx   = (xhigh - xlow) / mx;
        node = linspace(xlow-dx*mbc, xhigh+dx*mbc, mx + 2*mbc + 1);

        params.mx   = mx; params.node = node; params.dx   = dx;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % project the initial conditions.  note: meqn needs to be consistent with 
        % qinit, which is expected here.
        q   = zeros( mx, params.meqn, 1 );
        qex = zeros( mx, params.meqn, 1 );
        for n=1:mx
            q(n,:)    = qinit( node(n), node(n+1) ) / dx;
        end
        qex        = q;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        % take all the necessary time steps
        umax = params.u;
        dt   = dx / umax * CFL;
        t    = tstart;

        A = ConstructLinearLperiodic( params );
        params.A = A;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Main Time Integration loop                        %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nsteps = 0;
        q = reshape( q, mx*params.meqn, 1 );
        while( t < tfinal )

            if( t + dt > tfinal )
                dt = tfinal - t;
            else
                dt = dx / umax * CFL;
            end

            qn = q;
            q = rk_integrator( t, dt, qn );
            t = t + dt;

%           disp([['  dt = ', num2str(dt, '%2.5e' )]]);            
            nsteps = nsteps+1;

        end
        q = reshape( q, mx, params.meqn );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        % check error
        er(no) = norm( q - qex, 1 ) / norm( qex, 1 );

        % disp([['  er(', int2str(no),') = ', num2str(er(no),'%2.8e')]]);
        if( no > 1 )
            log2rat = log2( er(no-1) / er(no) );
            disp([[' log2( ratio of errors ) = ', num2str(log2rat, '%2.5f') ]] );
        end


    end
    % End of doing grid refinement %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp([['err =  ', num2str(er(1), '%2.10e' ) ]] );

    errp = er(2:end);
    errm = er(1:end-1);

    if( nrefine > 1 )
        disp([['log of the ratios = ']]);
        log2( errm./errp )
    end

    if( plt )

        % plot the solution
        xc = node( (1:mx) + mbc ) - dx/2;

        figure(1);
        clf;
        plot(xc, q(:,1), 'go');
        hold on;
        plot(xc, qex, 'b-' );
        hold off;

    end
