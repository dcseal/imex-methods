    clear all

    addpath('../lib/');

    format long;

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
            q(n,:)    = myquad( @qinit, node(n+mbc), node(n+1+mbc) ) / dx;
            aux(n,:)  = myquad( @auxinit, node(n+mbc), node(n+1+mbc) ) / dx;
            qex(n,:)  = myquad( @qexact, node(n+mbc), node(n+1+mbc) ) / dx;
        end
        q0         = q;
        params.aux = aux;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        % take all the necessary time steps
        umax = max( abs([ params.a, params.b] ) );
        dt   = dx / umax * CFL;
        t    = tstart;

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
            % q = sdc_integrator( t, dt, qn );
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
                    ';  \Delta t/\tau = ', num2str( dx / params.tau, '%2.3e' ), ...
                   '; error = ', num2str( er(no), '%2.5e' ), ...
                   '; last dt = ', num2str( dt, '%2.3e' ), ...
                   '; mx = ', num2str( mx, '%d' ), ...
                    '; nsteps = ', int2str(nsteps)]]);
        end
        
    end
    % End of doing grid refinement %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    errp = er(2:end);
    errm = er(1:end-1);

    if( plt )
        plot_results;
    end
