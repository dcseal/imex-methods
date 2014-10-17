    clear all

    set_params;

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
        qex  = zeros( length(tvec), meqn );
        qex(:,1)  = 1 ./ (1 + tvec );
        qex(:,2)  = tvec;

        q    = zeros( size(qex) );
        q(1,1) = 1.0;
        q(1,2) = tstart;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% Main Time Integration loop                        %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nsteps = 0;
        if( sdc )
            for n=1:mt
                % q(n+1,:)  = exp_integrator( tvec(n), dt, q(n,:) );
                q(n+1,:)  = ark_sdc_integrator( tvec(n), dt, q(n,:) );
                
                % Do not run the regular sdc_integrator as it wasn't designed to
                % handle implicit solves on delta when fI is non-linear.
                % q(n+1,:)  = sdc_integrator( tvec(n), dt, q(n,:) );

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





        % check error
        er(no) = norm( q - qex, 1 ) / norm( qex, 1 );
        if( no > 1 )
            log2rat = log2( er(no-1) / er(no) );
            disp([['; log2( ratio of errors ) = ', num2str(log2rat, '%2.5f'), ...
                   '; error = ', num2str( er(no), '%2.5e' ), ...
                   '; dt = ', num2str( dt, '%2.3e' ), '; mt = ', num2str(mt,'%3d')]] );
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
