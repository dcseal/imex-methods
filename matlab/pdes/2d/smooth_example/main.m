    clear all

    addpath('../lib/');

    % Read and print the parameters
    set_params;
    params


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Outer loop.  Each time this halves the grid size if nrefine > 1 %
    er = zeros(nrefine,1);
    for no=1:nrefine


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create new grid
        %
        mx   = mxstart*2^(no-1);  % number of grid cells (not counting boundary cells)
        my   = mystart*2^(no-1);  % number of grid cells (not counting boundary cells)

        params.mx = mx;
        params.my = my;


        % set up grid spacing and nodal points
        xhigh = params.xhigh;  xlow = params.xlow;  mbc = params.mbc;
        dx   = (xhigh - xlow) / mx;
        params.mx   = mx; params.dx   = dx;

        yhigh = params.yhigh;  ylow = params.ylow;  mbc = params.mbc;
        dy   = (yhigh - ylow) / my;
        params.my   = my; params.dy   = dy;

        nodex = linspace(xlow-dx*mbc, xhigh+dx*mbc, mx + 2*mbc + 1);
        nodey = linspace(ylow-dy*mbc, yhigh+dy*mbc, my + 2*mbc + 1);

%       node = zeros( (mx+1)*(my+1), 2 );
%       n = 1;
%       for i=1:(mx+1)
%       for j=1:(my+1)
%           node(n,1) = nodex(i);
%           node(n,2) = nodey(j);
%           n = n+1;
%       end
%       end
        nodex = linspace( xlow, xhigh, mx+1 );
        nodey = linspace( ylow, yhigh, my+1 );
        [xx yy] = meshgrid( nodex, nodey );

% TODO! 
params.node = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % project the initial conditions.  note: meqn needs to be consistent with 
        % qinit, which is expected here.
        tic();
        q   = zeros( mx, my, params.meqn, 1 );
        qex = zeros( mx, my, params.meqn, 1 );
        aux = zeros( mx, my, params.maux, 1 );
        for i=1:mx
        for j=1:my
            x = xx(j,i);  y = yy(j,i);
            q(i,j,:)    = myquad( @qinit, x, x+dx, y, y+dy );
            %aux(i,j,:)  = quad( @auxinit, node(n+mbc), node(n+1+mbc) ) / dx;
            %qex(i,j,:)  = qexact( node(n+mbc), node(n+1+mbc) ) / dx;
        end
        end
        q0         = q;
        qex = q0;
        params.aux = aux;

        %toc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        % take all the necessary time steps
        umax = params.u;
        dt   = dx / umax * CFL;
        t    = tstart;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Main Time Integration loop                        %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nsteps = 0;
        q = reshape( q, mx*my*params.meqn, 1 );
        while( t < tfinal )

            if( t + dt > tfinal )
                dt = tfinal - t;
            else
                dt = dx / umax * CFL;
            end

            qn = q;
            q = rk_integrator( t, dt, qn );
            t = t + dt;

            nsteps = nsteps+1;

        end
        q = reshape( q, mx, my, params.meqn );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        % check error
        er(no) = norm( q - qex, 1 ) / norm( qex, 1 );
        if( no > 1 )
            log2rat = log2( er(no-1) / er(no) );
            disp([[' log2( ratio of errors ) = ', num2str(log2rat, '%2.5f'), ...
                   ' er(', num2str(no,'%2d'),') = ', num2str(er(no),'%2.5e') ]] );
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
        plot_results;
    end
