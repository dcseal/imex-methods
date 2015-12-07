function ynp1 = semi_implicit_sdc_integrator( tn, dt, yn )
%Semi-implicit  SDC Integrator for the problem
%
%       y' = fI( y ) + fE( y ).
%
% Semi-implicit SDC updates the solution through
%
%   nu_{m+1} = nu_m + dt_m ( fE_{m  }^{k+1} - fE_{m  }^k )
%                     dt_m ( fI_{m+1}^{k+1} - fI_{m+1}^k ) + \int_m^{m+1} f^k dt.
%
% For each step of the method, the implicit equation that needs to be solve is
%
%       nu = nu_{exp} + dt_m fI( nu ),
% where
%
%       nu_exp = nu_m + dt_m ( fE_{m  }^{k+1} - fE_{m  }^k ) - fI_{m+1} + \int_m^{m+1} f^k dt.
%
% defines the "explicit" part of the integrator.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % coefficients necessary for performing integration
    [spts s] = get_quad_pts( params.sdc_order );
    coeffs = 0.5*dt*Res_Coeffs( spts );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a    = tn; b = tn+dt;
    tpts = 0.5* ( (a+b) + dt*spts );
    dt_vec = tpts(2:end) - tpts(1:(end-1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    k  = zeros( [s, length(yn) ] );
    kI = zeros( [s, length(yn) ] );
    kE = zeros( [s, length(yn) ] );

    % Intermediate storage for error and Q
    Q  = zeros( [s, length(yn) ] );
    Q(1,:) = yn';

    % 'Exact' value at time tn %
    kI(1,:) = fI  ( tn, yn )';
    kE(1,:) = fE  ( tn, yn )';

    % march q forward through each time step:
    for j=1:(s-1)

        qexp = Q(j,:) + dt_vec(j) * kE(j,:);

        % Take a Backward Euler step on each substep:
        Q  (j+1, :) = implicit_solve( qexp, tpts(j+1), dt_vec(j) );

        % Save the RHS (used for integrating the residual)
        kI (j+1, :) = fI  ( tpts(j+1), Q(j+1,:) )';
        kE (j+1, :) = fE  ( tpts(j+1), Q(j+1,:) )';

    end

    % Iterate over all the corrections
    for i=1:params.num_corrections

        % integrated residual
        N = coeffs * ( kI + kE );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Correct the solution
        for j=2:s

            % evaluate a new RHS part of the solution:
            tmp = k(j-1,:);

            kI(j-1, :) = fI  ( tpts(j-1), Q(j-1,:) )';
            kE(j-1, :) = fE  ( tpts(j-1), Q(j-1,:) )';

            qexp = Q(j-1,:) + dt_vec(j-1)*( kE(j-1,:) -tmp ) ...
                            - dt_vec(j-1)*  kI(j  ,:) + N(j-1,:);

            % Take a Backward Euler step on the new term:
            Q  (j, :) = implicit_solve( qexp, tpts(j), dt_vec(j-1) );

        end

        % Save the new right hand side (for evaluating integrals)
        kI(s, :) = fI  ( tpts(s), Q(s,:) )';
        kE(s, :) = fE  ( tpts(s), Q(s,:) )';
   
    end

    ynp1 = Q(s,:)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
