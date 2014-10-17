% Explicit right hand side function
function f = fE(t,y)

    global params
    mx   = params.mx;
    meqn = params.meqn;

    % source term gets treated implicitly.
    y = reshape( y, mx, meqn );
    f = ConstructL( y, t);
    f = reshape( f, mx*meqn, 1 );

end
