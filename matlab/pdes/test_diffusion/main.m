    clear all

    addpath('../lib/');

    set_params;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Outer loop.  Each time this halves the grid size if nrefine > 1 %
    er = zeros(nrefine,1);
    for no=1:nrefine


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create new grid
        %
        mx   = mxstart*2^(no-1);  % number of grid cells (not counting boundary cells)
        params.mx = mx;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% TODO - remove me %%%%
%       mx = mxmax;
%       params.mx = mxmax;
%       CFL = 0.5*CFL;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        aux = zeros( mx, params.maux, 1 );
        for n=(1:mx)
            q(n,:)    = qinit( node(n+mbc), node(n+1+mbc) ) / dx;
            aux(n,:)  = quad( @auxinit, node(n+mbc), node(n+1+mbc) ) / dx;
            qex(n,:)  = qexact( node(n+mbc), node(n+1+mbc) ) / dx;
        end
        q0         = q;
        params.aux = aux;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % figure out length of desired time step
        umax = max( params.u, 1.0 );
        t    = tstart;
        mt = ceil( (tfinal-tstart) * CFL / (umax*dx) );
        dt = (tfinal -  tstart ) / mt;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Main Time Integration loop                        %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nsteps = 0;
        q = reshape( q, mx*params.meqn, 1 );
        for n=1:mt

            qn = q;
            % q = rk_integrator( t, dt, qn );
            q = sdc_integrator( t, dt, qn );
            t = t + dt;

            nsteps = nsteps+1;

        end
        q = reshape( q, mx, params.meqn );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        % check error
        er(no) = norm( q - qex, 1 ) / norm( qex, 1 );
        if( no > 1 )
            log2rat = log2( er(no-1) / er(no) );
            disp([['; log2( ratio of errors ) = ', num2str(log2rat, '%2.5f'), ...
                   '; error = ', num2str( er(no), '%2.5e' ), ...
                   '; Delta t / tau = ', num2str( dt / params.eps, '%2.3e' ), ...
                   '; mx = ', num2str( mx, '%d' ), ...
                   '; mt = ', num2str( mt, '%d' )]] );
        end
        
    end
    % End of doing grid refinement %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if( plt )
        plot_results;
    end
