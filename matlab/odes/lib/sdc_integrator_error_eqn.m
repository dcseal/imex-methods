function ynp1 = sdc_integrator( tn, dt, yn )
% SDC Integrator based on forward and backward Euler steps %.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [spts s] = get_weights( params.sdc_order );

    % coefficients necessary for performing integration
    coeffs = 0.5*dt*Res_Coeffs( spts );

    a    = tn; b = tn+dt;
    tpts = 0.5* ( (a+b) + dt*spts );
    dt_vec = tpts(2:end) - tpts(1:(end-1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    kI = zeros( [s, length(yn) ] );
    kE = zeros( [s, length(yn) ] );

    % Intermediate storage for error and Q
    de = zeros( [s, length(yn) ] );
    Q  = zeros( [s, length(yn) ] );
    Q(1,:) = yn';

    % 'Exact' value at time tn %
    kI(1,:) = fI( tn, yn )';
    kE(1,:) = fE( tn, yn )';


    % march q forward through each time step:
    for j=1:(s-1)
        qexp = Q(j,:) + dt_vec(j) * kE(j,:);
        Q  (j+1, :) = implicit_solve( qexp, tpts(j+1), dt_vec(j) );
        kE (j+1, :) = fE( tpts(j+1), Q(j+1,:) )';
        kI (j+1, :) = fI( tpts(j+1), Q(j+1,:) )';
    end

    % Iterate over all the corrections
    for i=1:params.num_corrections

        % evaluate the new right hand side in order to form the residual
        if( i > 1 )
            for j=2:s
                kE (j, :) = fE( tpts(j), Q(j, :) )';
                kI (j, :) = fI( tpts(j), Q(j, :) )';
            end
        end

        % integrated residual
        N = coeffs * ( kE + kI );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % march the error forward through each sub-interval:

        % Note: first step is special, since de(1) = 0;  We can cut out two
        % function evaulations here 
        % explicit step
        j = 2;
        de(j,:) = - ( Q(j,:)-Q(j-1,:) ) + N(j-1,:);

        % implicit step:
        de(j,:) = implicit_solve( de(j,:), tpts(j), dt_vec(j-1) );

        % remaining sub-intervals
        for j=3:s

            % TODO - for some reason this doesn't work well with the stability 
            % plot.  when i replace kE with "fE(0, Q(j-1,:) )" , I get 
            % reasonable plots

            % explicit step
            de(j,:) = de(j-1,:) - ( Q(j,:)-Q(j-1,:) ) ...
               + dt_vec(j-1) * ( ...
                fE( tpts(j-1), Q(j-1,:) + de(j-1,:) ) - kE( j-1, : ) ) ...
               + N(j-1,:);

            % implicit step:
            de(j,:) = implicit_solve( de(j,:), tpts(j), dt_vec(j-1) );

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
% This function determines the quadrature points used for integration of the
% residual.  The default behaviour is to use uniform time points.

    if( N == 4 )

        % These weights and points integrate a polynomial of degree 5 exactly,
        % and therefore are 6th order accurate.
%       spts = [ -1.00000000000000000000e+00, -4.47213595499957927704e-01, ...
%                4.47213595499957927704e-01, 1.00000000000000000000e+00];
%       spts = [-1, -1/sqrt(3), 1/sqrt(3), 1];
%       spts = [-1 -0.5 0.5 1];
        spts = linspace(-1,1,4);

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
