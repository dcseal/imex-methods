function ynp1 = backward_taylor_imex( tn, dt, yn )
%EXPERIMENTAL Backward Taylor IMEX SOLVER
%
%       y' = A y + fE( y ).
%
% This scheme takes a backward Taylor expansion of y and replaces each
% explicit piece of the backward Taylor series with a forward Taylor
% expansion.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Matrix for implicit portion of the solver
    A = params.A;
    I = eye( length(yn) );

    fn = fE( tn, yn );
    Jn = jacob( tn, yn );

    % Matrix to invert
    M    = I - dt*A + dt^2/2*( A^2 - Jn*A );

    % explicit portion of step
    yexp = yn + dt*fn + dt^2/2*( Jn*fn - A*fn );

    % Full time step
    ynp1 = M \ yexp;

end
