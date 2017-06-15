%%% Main program.  See set_params for user defined options. %%%

    clear all;

    addpath('../lib/');

    set_params;

    if( isfield( params, 'fname' ) )
        fids = fopen(params.fname, 'w');
    end

    meqn = params.meqn;

    % make the output directory
%   cmd       = ['mkdir ', outputdir];
%   [status result] = system(cmd);
%   if( status )
%       disp([ 'Warning: unable to make directory: ', outputdir] );
%       disp(result)
%   end

    % print a helpful output file describing this run
%   fname = [outputdir, 'qhelp.dat'];
%   fid   = fopen(fname, 'w');
%   fprintf(fid,'%d\n', nrefine);
%   fprintf(fid,'%d\n', meqn);
%   fprintf(fid,'%2.15e\n', tstart);
%   fprintf(fid,'%2.15e\n', tfinal);
%   fprintf(fid,'%2.15e\n', params.eps);
%   fclose(fid);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Outer loop.  Each time this halves the grid size if nrefine > 1 %
    er = zeros(nrefine,1);
    for no=1:nrefine

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mt   = mstart * (2^(no-1));
        dt   = (tfinal - tstart ) / mt;
        t    = tstart;
        tvec = linspace( tstart, tfinal, mt+1 );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initial Conditions (note: there is no exact solution for this problem
        q        = zeros( mt+1, meqn );
        q( 1, 1 ) = 2;
        q( 1, 2 ) = 2./3.;
        q( 1, 2 ) = -0.666666654321;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        opts = odeset('RelTol', 1e-13, 'AbsTol',1e-14);
%       [tvals,qex] = ode23(@vdp1,tvec,[q(1,1); q(1,2)], opts);
        [tvals,qex] = ode15s(@vdp1,tvec,[q(1,1); q(1,2)], opts);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Main Time Integration loop                        %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic();
        nsteps = 0;
        for n=1:mt

            qn = q(n,:);
            qn = advance( tvec(n), dt, qn );
            q( n+1, :)  = qn;
            nsteps = nsteps+1;

        end
        run_time = toc();
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if( isfield( params, 'fname' ) )
            er1 = norm( q - qex, 1 ) / norm( qex, 1 );
            er2 = norm( q - qex, 2 ) / norm( qex, 2 );
            eri = max( abs( q - qex ) ) / max( abs( qex ) );
            fprintf(fids, '%i %2.15e %2.15e %2.15e\n', mt, er1, er2, eri );
            fprintf(1, '%i %2.15e %2.15e %2.15e\n', mt, er1, er2, eri );
        end

        % print output to file
%       tic()
%       fname = [outputdir, 'q', num2str(no,'%04d'), '.dat'];
%       fid   = fopen(fname, 'w');
%       fprintf(fid,'%d\n', mt);
%       for n=1:(mt)
%           for m=1:meqn
%               fprintf(fid, '%2.15e ', q(n,m) );
%           end
%           fprintf(fid, '\n');
%       end 
%       fclose(fid);
%       disp([['It took ', num2str( toc(), '%2.3e' ), '  seconds to print to file  ', ...
%              'It took ', num2str( run_time, '%2.3e'  ), '  seconds run the integration loop']]);

    end
    % End of doing grid refinement %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if( isfield( params, 'fname' ) )
        fclose( fids );
    end


    if( plt )
        figure(1);
        clf;

        q1 = q(:,1);
        q2 = q(:,2);

        plot( tvec, q1, 'go' );
        hold on;
        plot( tvec, q2, 'ro' );
        plot( tvec, qex );
        hold off;
        
        l1 = legend('y1', 'y2', 'Exact' );
        set( l1, 'FontSize', 16 );

        ttl = ['Van der Pol.  tfinal = ', num2str( tfinal ), '  \epsilon = ', ...
             num2str( params.eps, '%2.1e') ];
        ttl = 'Van der Pol.';

        t = title( ttl );
        set( t, 'FontSize', 16 );

        axis([0 tfinal -5 3]);
    end
