%%%%%%%%%% Kaps main %%%%%%%%%%%%%

    clear all

    addpath('../lib/');

    set_params;

    meqn = params.meqn;

    % make the output directory
    cmd       = ['mkdir ', outputdir];
    [status result] = system(cmd);
    if( status )
        disp([ 'Warning: unable to make directory: ', outputdir] );
        disp(result)
    end

    % print a helpful output file describing this run
    fname = [outputdir, 'qhelp.dat'];
    fid   = fopen(fname, 'w');
    fprintf(fid,'%d\n',     nrefine);
    fprintf(fid,'%d\n',     meqn);
    fprintf(fid,'%2.15e\n', tstart);
    fprintf(fid,'%2.15e\n', tfinal);
    fprintf(fid,'%2.15e\n', params.tau);
    fclose(fid);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Outer loop.  Each time this halves the grid size if nrefine > 1 %
    er = zeros(nrefine,1);
    for no=1:nrefine


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mt = mstart*2^(no-1);
        dt   = (tfinal - tstart ) / mt;
        t    = tstart;
        tvec = linspace( tstart, tfinal, mt+1 );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        qex = 1/(1+params.tau^2)*( cos(tvec) + params.tau * sin(tvec) ...
                - exp( -tvec / params.tau ) )';
        q  = zeros( size(qex) );
        q0 = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Main Time Integration loop                        %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nsteps = 0;
        if( sdc )
            for n=1:mt
                q(n+1,:)  = ark_sdc_integrator( tvec(n), dt, q(n,:) );
                %q(n+1,:)  = sdc_integrator( tvec(n), dt, q(n,:) );
                nsteps = nsteps+1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            for n=1:mt
                q(n+1,:)  = rk_integrator( tvec(n), dt, q(n,:) );
                nsteps = nsteps+1;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % print output to file
        fname = [outputdir, 'q', num2str(no,'%04d'), '.dat'];
        fid   = fopen(fname, 'w');
        fprintf(fid,'%d\n', mt);
        for n=1:(mt)
            for m=1:meqn
                fprintf(fid, '%2.15e ', q(n,m) );
            end
            fprintf(fid, '\n');
        end 
        fclose(fid);



        % check error
        er(no) = norm( q - qex, 1 ) / norm( qex, 1 );
        if( no > 1 )
            log2rat = log2( er(no-1) / er(no) );
            disp([['; log2( ratio of errors ) = ', num2str(log2rat, '%2.5f'), ...
                   '; error = ', num2str( er(no), '%2.5e' ), ...
                   '; dt = ', num2str( dt, '%2.3e' ), '; \Delta t / tau = ', ...
                    num2str(dt/params.tau, '%2.3e')]] );
        end
        
    end
    % End of doing grid refinement %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp([['err =  ', num2str(er(1), '%2.10e' ) ]] );

    errp = er(2:end);
    errm = er(1:end-1);

    if( plt )
        plot_results;
    else
        exit;
    end
