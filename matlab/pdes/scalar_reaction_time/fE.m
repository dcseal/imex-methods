% Explicit right hand side function
function f = fE(t,y)

    global params
    mx   = params.mx;
    meqn = params.meqn;
    aux  = params.aux;

    % advection term plus source term.
    y = reshape( y, mx, meqn );
    f = ConstructL( y, t) + cos(t)*aux/params.tau;
    f = reshape( f, mx*meqn, 1 );

end
