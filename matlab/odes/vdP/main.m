%%% Main program.  See set_params for user defined options. %%%

    clear all;

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
    fprintf(fid,'%d\n', nrefine);
    fprintf(fid,'%d\n', meqn);
    fprintf(fid,'%2.15e\n', tstart);
    fprintf(fid,'%2.15e\n', tfinal);
    fprintf(fid,'%2.15e\n', params.eps);
    fclose(fid);



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
        q( 1, 2 ) = -0.6666654321121172;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Main Time Integration loop                        %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic();
        nsteps = 0;
        for n=1:mt

            qn = q(n,:);
            % qn = rk_integrator( tvec(n), dt, qn );
            % qn = ark_sdc_integrator( tvec(n), dt, qn );
            qn = semi_implicit_sdc_integrator( tvec(n), dt, qn );
            q( n+1, :)  = qn;
            nsteps = nsteps+1;

        end
        run_time = toc();
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % print output to file
        tic()
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
        disp([['It took ', num2str( toc(), '%2.3e' ), '  seconds to print to file  ', ...
               'It took ', num2str( run_time, '%2.3e'  ), '  seconds run the integration loop']]);

    end
    % End of doing grid refinement %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if( plt )
        plot_results;
    else
        exit;
    end
