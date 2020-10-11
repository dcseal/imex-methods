function ynp1 = backward_taylor_imex( tn, dt, yn )
%EXPERIMENTAL Backward Taylor IMEX SOLVER
%
%       y' = A y + f(t)
%
% This scheme takes a backward Taylor expansion of y and replaces each
% explicit piece of the backward Taylor series with a forward Taylor
% expansion.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Matrix for implicit portion of the solver
    A = -1/params.tau;
    I = eye( length(yn) );

    fn  = fE( tn, yn );
    dfn = -sin(t)/params.tau;

    % Matrix to invert
    M    = I - dt*A + dt^2/2*( A^2 );

    % explicit portion of step
    yexp = yn + dt*fn + dt^2/2*( dfn - fn );

    % Full time step
    ynp1 = M \ yexp;

end
