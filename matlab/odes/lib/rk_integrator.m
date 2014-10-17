% RK Integrator.  This routine takes a single time step using the butcher
% tableau specified by the file being run.
function ynp1 = rk_integrator( tn, dt, yn, de )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Grab the Butcher tableau %%%%%%%%%

    coef = params.coeffs;
    if(     strcmp( coef, 'coeffs_imex' ) )
        run coeffs_imex
    elseif( strcmp( coef, 'coeffs_ark32' ) )
        run coeffs_ark32
    elseif( strcmp( coef, 'coeffs_ark43' ) )
        run coeffs_ark43
    elseif( strcmp( coef, 'RK2' ) )
        run coeffs_classicalRK2;
    elseif( strcmp( coef, 'RK3' ) )
        run coeffs_classicalRK3;
    elseif( strcmp( coef, 'coeffs_ark2ars' ) )
        run coeffs_ark2ars;
    else
        run coeffs_classicalRK4;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % storage for all the intermediate time points
    kI = zeros( [s, length(yn) ] );
    kE = zeros( [s, length(yn) ] );

    for i=1:s

        % `explicit' part of the time step
        yi = yn;
        for j=1:i-1
            yi = yi + dt * ( AE(i,j)*kE(j,:) + AI(i,j)*kI(j,:) );
        end

        % implicit solve:
        if( i > 1 )

            % parameters needed for the implicit solve:
            % yi = (yexp)_i + as * fI(t, y_i) %
            global yexp as ts;
            ts   = tn+dt*c(i);
            as   = dt*AI(i,i);
            yexp = yi;
            yi = implicit_solve( yi, ts, as );

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Compute intermediate stage values    %%%
        kE(i,:) = fE( tn + c(i)*dt, yi )';
        kI(i,:) = fI( tn + c(i)*dt, yi )';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update the solution
    ynp1 = yn;
    for n=1:s
        ynp1 = ynp1 + dt*( bE(n)*kE(n,:) + bI(n)*kI(n,:) );
    end

    ynp1 = ynp1';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end
