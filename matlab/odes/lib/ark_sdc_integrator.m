function ynp1 = ark_sdc_integrator( tn, dt, yn )
% SDC Integrator based on forward and backward Euler steps %.


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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [spts s] = get_weights( params.sdc_order );

    % coefficients necessary for performing integration
    coeffs = 0.5*dt*Res_Coeffs( spts );

    a    = tn; b = tn+dt;
    tpts = 0.5* ( (a+b) + dt*spts );
    dt_vec = tpts(2:end) - tpts(1:(end-1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    kE = zeros( [s, length(yn) ] );

    % Intermediate storage for error and Q
    de = zeros( [s, length(yn) ] );
    Q  = zeros( [s, length(yn) ] );
    Q(1,:) = yn';

    % 'Exact' value at time tn %
    kE(1,:) = fE( tn, yn )';
    kI(1,:) = fI( tn, yn )';

    % march q forward through each time step:
    for j=1:(s-1)

        h = dt_vec(j);
        local_tpts = tpts(j) + c * h;

        Q(j+1,:)  = rk_integrator( tpts(j), h, Q(j,:) );
        kE(j+1,:) = fE( tpts(j+1), Q(j+1,:) )';
        kI(j+1,:) = fI( tpts(j+1), Q(j+1,:) )';

    end

    % Iterate over all the corrections
    num_corrections = params.num_corrections;
    for i=1:num_corrections
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % march the error forward through each sub-interval:

        if( i > 1 )
            for j=2:s
                kE(j,:) = fE( tpts(j), Q(j,:) );
                kI(j,:) = fI( tpts(j), Q(j,:) );
            end
        end

        for j=1:(s-1)

            h = dt_vec(j);
            local_tpts = tpts(j) + c * h;
            local_spts = spts(j) + (spts(j+1)-spts(j))*c;
            [EvalCoeffs Integrate_Coeffs] = New_Res_Coeffs( spts, local_spts );

            FEinterp     = EvalCoeffs * kE;
            FIinterp     = EvalCoeffs * kI;
            Qinterp      = EvalCoeffs * Q;
            Nintegrated  = 0.5*dt*Integrate_Coeffs * ( kE + kI );

            de(j+1,:) = rk_integrator_extra(                ...
                            tpts(j), h, de(j,:),            ...
                        Qinterp, Q(j,:), FEinterp, kE(j,:), ...
                        FIinterp, kI(j,:), Nintegrated );

            %%%%%%%%%% Classical Explicit RK4 Time Step %%%%%%%%%%%%%%
%           qs = Q(j,:);

%           d1 = de(j,:);
%           k1 = fE( local_tpts(1), qs + d1 ) - fE(local_tpts(1), qs);

%           d2 = d1-(Qinterp(1,:)-Q(j,:)) + Nintegrated(1,:) + 0.5*h*k1;
%           k2 = fE( local_tpts(2), d2 + Qinterp(1,:) ) - fE(local_tpts(2), Qinterp(1,:) );

%           d3 = d1-(Qinterp(2,:)-Q(j,:)) + Nintegrated(2,:) + 0.5*h*k2;
%           k3 = fE( local_tpts(3), d3 + Qinterp(2,:) ) - fE(local_tpts(3), Qinterp(2,:) );

%           d4 = d1-(Qinterp(3,:) - Q(j,:)) + Nintegrated(3,:) + 1.0*h*k3;
%           k4 = fE( local_tpts(4), d4 + Qinterp(3,:) ) - fE( local_tpts(4), Qinterp(3,:) );

%           de(j+1,:) = de(j,:) - (Q(j+1,:)-Q(j,:)) + Nintegrated(3,:) ...
%               + h * ( bE(1)*k1 + bE(2)*k2 + bE(3)*k3 + bE(4)*k4 );
 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % add the error into Q:
        for j=2:s
            Q(j,:) = de(j,:) + Q(j,:);
        end

    end
    ynp1 = Q(s,:)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function [spts s] = get_weights(N)

    if( N == 4 )

        % These weights and points integrate a polynomial of degree 5 exactly,
        % and therefore are 6th order accurate.
%       spts = [ -1.00000000000000000000e+00, -4.47213595499957927704e-01, ...
%                4.47213595499957927704e-01, 1.00000000000000000000e+00];
        spts = linspace(-1,1,4);
%       spts = [-1, -1/sqrt(3), 1/sqrt(3), 1];

        wgts = [1.66666666666666657415e-01, 8.33333333333333370341e-01, ...
                8.33333333333333370341e-01, 1.66666666666666657415e-01];

    elseif( N == 3 )

        % These weights and points integrate a polynomial of degree 3 exactly,
        % and therefore are 4th order accurate.
        spts = [-1.0, 0.0, 1.0];

        wgts = [ 3.33333333333333314830e-01, ...
                 1.33333333333333325932e+00, ...
                 3.33333333333333314830e-01];

    elseif( N == 2 )
    
        % 2nd Order SDC
        spts = [-1 1];
        wgts = [1 1];

    else
        spts = linspace(-1,1,N);
    end

    s = length( spts );

end
