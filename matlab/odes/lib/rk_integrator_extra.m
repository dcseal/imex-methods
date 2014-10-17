% RK Integrator.  This routine takes a single time step using the butcher
% tableau specified by the file being run.
function dnp1 = rk_integrator_extra( tn, dt, dn, ...
    Qinterp, qs, FEinterp, FEstart, ...
    FIinterp, FIstart, Nintegrated )

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
    kI = zeros( [s, length(dn) ] );
    kE = zeros( [s, length(dn) ] );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i = 1;
    kE(i,:) = fE( tn + c(i)*dt, dn + qs ) - FEstart;
    kI(i,:) = fI( tn + c(i)*dt, dn + qs ) - FIstart;
    for i=2:(s-1)

        % `explicit' part of the time step
        de_i = dn + Nintegrated(i-1,:) - (Qinterp(i-1,:)-qs);
        for j=1:i-1
            de_i = de_i + dt * ( AE(i,j)*( kE(j,:) ) + AI(i,j)*kI(j,:) );
        end

        % parameters needed for the implicit solve:
        % de_i = (yexp)_i + as * fI(t, y_i) %
        global yexp as ts;
        ts   = tn+dt*c(i);
        as   = dt*AI(i,i);
        yexp = de_i;
        de_i = implicit_solve_for_delta( de_i, Qinterp(i-1,:), ts, as );


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Compute intermediate stage values    %%%
        kE(i,:) = fE( tn + c(i)*dt, de_i+Qinterp(i-1,:) ) - FEinterp(i-1,:);
        kI(i,:) = fI( tn + c(i)*dt, de_i+Qinterp(i-1,:) ) - FIinterp(i-1,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    i=s;

    % `explicit' part of the time step
    de_i = dn + Nintegrated(i-1,:) - (Qinterp(i-1,:)-qs);
    for j=1:i-1
        de_i = de_i + dt * ( AE(i,j)*( kE(j,:) ) + AI(i,j)*kI(j,:) );
    end

    % parameters needed for the implicit solve:
    % de_i = (yexp)_i + as * fI(t, y_i) %
    global yexp as ts;
    ts   = tn+dt*c(i);
    as   = dt*AI(i,i);
    yexp = de_i;
    de_i = implicit_solve_for_delta( de_i, Qinterp(i-1,:), ts, as );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute intermediate stage values    %%%
    kE(i,:) = fE( tn + c(i)*dt, de_i+Qinterp(s-1,:) ) - FEinterp(s-1,:);
    kI(i,:) = fI( tn + c(i)*dt, de_i+Qinterp(i-1,:) ) - FIinterp(i-1,:);
    % kI(i,:) = fI( tn + c(i)*dt, de_i );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update the solution
    dnp1 = dn + Nintegrated(s-1,:) - ( Qinterp(s-1,:) - qs );
    for n=1:s
        dnp1 = dnp1 + dt*( bE(n)*kE(n,:) + bI(n)*kI(n,:) );
    end

    dnp1 = dnp1';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
