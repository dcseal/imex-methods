% RK Integrator.  This routine takes a single time step using the butcher
% tableau specified by the file being run.
function ynp1 = rk_integrator( tn, dt, yn )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Grab the Butcher tableau %%%%%%%%%

    % option 1:
    run coeffs_imex

    % option 2:
    % run coeffs_classicalRK4

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % storage for all the intermediate time points
    kI = zeros( [s, size(yn) ] );
    kE = zeros( [s, size(yn) ] );



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global params
    mx   = params.mx;
    meqn = params.meqn;
%   r    = params.r;
%   tau  = params.tau;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    for i=1:s

        % `explicit' part of the time step
        yi = yn;
        for j=1:i-1
            yi = yi + dt * ( AE(i,j)*kE(j,:)' + AI(i,j)*kI(j,:)' );
        end

        % parameters needed for the implicit solve:
        % yi = (yexp)_i + as * fI(t, y_i) %
        global yexp as ts;
        ts   = tn+dt*c(i);
        as   = dt*AI(i,i);
        yexp = yi;



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % implicit solve for the current value of y
%       [yi, fval, exitflag] = fsolve( @fImplicitFunc, yi );
%       if( exitflag > 1 || exitflag < 1 )
%           disp('fsolve did not converge.  garbage is being computed');
%           disp([['exitflag = ', int2str(exitflag)]]);
%       end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % implicit solve here:
        yi = ( eye( mx ) - params.uimp * as*params.A ) \  yexp;
%       yi = implicit_solve( yi, ts, as );


        % First problem (meqn = 1):
%       yi = ( yexp - as/tau*aux ) / (1+as/tau);




        % Second problem (meqn = 2, stiff source):
%       yi = reshape(yi, mx, meqn );
%       yi(:,2) = ( yi(:,2) + as/tau*r*yi(:,1) ) / ...
%                 ( 1+as/tau );

%       yi = reshape( yi, mx*meqn, 1 );




        % Third problem (meqn = 1, scale difference):
%       yi = reshape(yi, mx, meqn );
%       A = [1+as/tau, -as/tau; -as/tau, 1+as/tau];
%       [yi(:,1) yi(:,2)] = A \ [yi(:,1) yi(:,2)];
%       yi = reshape( yi, mx*meqn, 1 );




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Compute intermediate stage values    %%%

        % evaluate the right hand side at this time value
        kE(i,:) = fE( tn + c(i)*dt, yi )';

        % TODO - there should be a way to get kI from the implicit solve above!
        kI(i,:) = fI( tn + c(i)*dt, yi )';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update the solution
    ynp1 = yn';
    for n=1:s
        ynp1 = ynp1 + dt*( bE(n)*kE(n,:) + bI(n)*kI(n,:) );
    end

    ynp1 = ynp1';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end
